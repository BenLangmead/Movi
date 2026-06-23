# SSHash-compatible minimal perfect hash (MPH) for Movi

Goal: expose, from a single Movi index, a k-mer → integer-id dictionary whose
interface and semantics match SSHash's MPH, so Movi can be a drop-in for code that
uses SSHash's `lookup` / `access` (e.g. attaching per-k-mer satellite data).

This builds on the existing `build-kmerbv` / `--kmer --kmer-bv` work (per-k
bitvector `B_k` marking the first BWT row of each distinct real k-mer interval,
with `rank(lb)` giving a dense id).

## SSHash's MPH interface (jermp/sshash `include/dictionary.hpp`, `util.hpp`)

```cpp
// forward: k-mer -> id (+ metadata)
lookup_result lookup(char const* string_kmer, bool check_reverse_complement = true) const;
lookup_result lookup(Kmer  uint_kmer,        bool check_reverse_complement = true) const;
bool          is_member(...) const;

// inverse: id -> k-mer
void          access(uint64_t kmer_id, char* string_kmer) const;

// cardinality / params
uint64_t num_kmers()  const;     // size of the id space; ids are in [0, num_kmers)
uint64_t k()          const;
uint64_t m()          const;     // minimizer length (N/A for Movi)
bool     canonical()  const;

// iteration over (id, kmer)
struct iterator { bool has_next(); std::pair<uint64_t, Kmer> next(); };
iterator begin() const;          // and at_kmer_id(id)
```

`lookup_result` (the fields that are meaningful for an MPH):

```cpp
struct lookup_result {
    uint64_t kmer_id;          // absolute id in [0, num_kmers); invalid_uint64 if absent
    uint64_t kmer_id_in_string;// SPSS-relative id      -- SSHash-specific, no Movi analog
    uint64_t kmer_offset;      // SPSS offset           -- SSHash-specific
    int64_t  kmer_orientation; // forward_orientation (0) or backward_orientation (1)
    uint64_t string_id, string_begin, string_end;  // SPSS contig -- SSHash-specific
    bool     minimizer_found;  // SSHash-specific
};
```

Key semantics (with `--canonical`, which is what we compare against): **one id per
canonical k-mer**, ids dense in `[0, num_kmers)`. `lookup(x)` and `lookup(rc(x))`
return the **same** `kmer_id` and opposite `kmer_orientation`. `access(id)` returns
the canonical (forward-orientation) spelling.

## The semantic gap vs. Movi's current kmer-bv

Movi's index is the doubled text `T = ref $ rc(ref)`. Searching the literal pattern
`x` over `T` yields `occ_T(x) = occ_ref(x) + occ_ref(rc(x))` — the **canonical
count** (already correct, matches KMC). But the *interval position* differs for `x`
vs `rc(x)`: they are two distinct lexicographic patterns with two distinct `lb`
rows, hence two distinct `rank(lb)` ids. So today:

- **count** is canonical ✓
- **id** is **per-orientation**: `lookup(x) != lookup(rc(x))`, and the id space has
  size ≈ `2 · (#canonical k-mers)` (minus palindromes), **not** `num_kmers`.

For SSHash compatibility we need one id per *canonical* k-mer in `[0, num_kmers)`.

## Design

### 1. Canonical id space (`build-kmerbv --canonical`)

In `build_kmerbv`, the DFS reaches every distinct real k-mer once as a depth-`k`
node. Track the k-mer **string** down the DFS (prepend each branch base), and at a
depth-`k` node set `B_k[lb] = 1` **only if the k-mer is canonical** (`s == min(s,
rc(s))`). Then:

- `rank_1(B_k)` counts only canonical k-mers ⇒ ids dense in `[0, num_canonical)`,
  matching SSHash's `num_kmers`.
- non-canonical orientations are simply not marked.

Keep the existing per-orientation mode (no `--canonical`) for back-compat and for
the counting path (which is orientation-agnostic anyway). Record the mode in a
1-byte sidecar (`kmerbv.<k>.meta`) so the query side knows whether ids are
canonical.

### 2. `lookup` (forward)

```
lookup(x):
  c = canonical(x) = lexicographically-min(x, rc(x))   # only when canonical mode
  iv = backward_search(c)                               # single k-mer interval
  if iv empty: return { kmer_id = invalid, found=false }
  lb = all_p[iv.run_start] + iv.offset_start
  return { kmer_id = rank_1(B_k, lb),
           kmer_orientation = (x == c ? forward : backward),
           occ_count = iv.count(rlbwt) }               # bonus over SSHash
```

Maps onto `lookup_result.{kmer_id, kmer_orientation}`. The SPSS-specific fields
(`string_*`, `kmer_id_in_string`, `minimizer_found`) have no Movi analog and are
left `invalid` — documented.

### 3. `access` (inverse, id → k-mer)

Add `select_support_mcl` over `B_k`. Then:

```
access(id):
  lb = select_1(B_k, id + 1)             # first row of the id-th canonical k-mer
  reconstruct k chars of the suffix at row lb by forward-stepping the move
  structure k times (Movi's forward/inverse-LF navigation), emitting the
  canonical spelling.
```

This is the one genuinely new traversal; it needs Movi forward navigation from a
BWT row. `num_bits` cost: one `select` structure per k (only built/loaded on
demand for `access`).

### 4. Cardinality / params

- `num_kmers(k)` = `rank_1(B_k, B_k.size())` (already computed at build; store it).
- `k()` from the query `-k`; `canonical()` from the sidecar; `m()` is N/A (Movi has
  no minimizer) → report 0 / "n/a".

### 5. CLI / output

- Build: `movi build-kmerbv --kmer-lengths k1,k2,... [--canonical]`.
- Lookup: `movi query --kmer --kmer-bv -k K` → per k-mer emit `id[:orientation]`
  (and `:count`), already wired through `mq.add_kmer(pos, found, kmer_id, occ)`;
  add the orientation field and a header line documenting the id space.
- Access: `movi access-kmer --kmer-bv -k K --ids <file>` → one k-mer per id.
- A `--print-info` that prints `num_kmers`, `k`, `canonical` like SSHash's
  `print_info()`.

## Validation

- **Bijection**: for every canonical k-mer enumerated from the index, `access(lookup(x).kmer_id) == x`, and ids cover `[0, num_kmers)` exactly once (extend the existing k↔id bijection check to canonical mode at k=15/21/31).
- **Cross-tool**: against an SSHash index built `--canonical` on the same input,
  `num_kmers` matches, and (since ids are tool-specific orderings) the *set*
  `{ (canonical-kmer) }` and per-k-mer membership/count agree. Id *values* need not
  match SSHash's (different MPHF), but must be a valid dense bijection — that is the
  compatibility contract (same as any two distinct MPHFs).

## Status / plan

- [x] `build-kmerbv --canonical`: 2-bit packed k-mer carried down the DFS, mark `B_k` only when `packed <= rc(packed)`; `kmerbv.<k>.meta` sidecar (canonical flag + num_kmers). **Validated: ecoli10 k=31 num_kmers=15,209,176 == KMC unique exactly; non-canonical = 2× (no odd-k palindromes).** (commit `a500502`)
- [x] **sd_vector is the default id rep; dense is opt-in** (`build-kmerbv --id-bv-dense`). Default build writes only `kmerbv.<k>.sd`; load defaults to sd, `MOVI_ID_BV=dense` forces dense. (validated: ecoli10 default build writes only `.sd`, query loads it.)
- [x] **`lookup` canonical ids — correctness DONE & VALIDATED.** In canonical mode `--kmer --kmer-bv` now returns `id = rank(min(lb_fw, lb_rc))`. **Bijection test (ecoli10 k=31, num_kmers=15,209,176): max id 15,209,173 < num_kmers, 0 out-of-range, and the forward-reads vs reverse-complement-reads id multisets are IDENTICAL ⇒ `lookup(x)==lookup(rc(x))`.** This is the SSHash-compatible canonical id space, queryable.
- [x] **orientation output** — `--kmer --kmer-bv` now emits `pos:count:id:f|r` (`add_kmer` 4th field), matching SSHash's `lookup_result.kmer_orientation`; set from `lb_fw <= lb_rc`. Validated (`:f`/`:r` present, bijection still holds).
- [x] **`lookup` speed — lazy-rc (DONE).** Canonicality is decided from the strings (`x <= rc(x)`, == `lb_fw <= lb_rc`), so canonical k-mers use `lb_fw` directly with NO rc work; only non-canonical k-mers compute the rc interval, via the normal **ftab-assisted** `initialize_backward_search`+`backward_search` on the rc string. **Result: 0.75 M/s vs 0.48 M/s always-rc = +56% (1.56×), byte-identical ids** (ecoli100 k31, 60 M kmers). Canonicalization now costs ~25% over the non-canonical 1.0 M/s (was 50%). Now the default (no env gate). The bidirectional keep-going walk below would push further (toward the count-bv ~3–5 M/s) but is a larger rewrite; lazy-rc captures most of the win cheaply.
- [x] **`lookup` keep-going id walk — DONE & DEFAULT.** `query_kmers_id_bv` mirrors the count-bv positive-skip walk but emits the canonical MPHF id per k-mer using only `B_k`+`all_p` (no count structure): for a present k-mer the walk's subset interval gives `id = rank_B(lb_row+1)-1` (x's mark is the only `B_k` 1 in its group); non-canonical k-mers use the lazy ftab rc search. **Result: 1.61 M/s vs 0.73 M/s per-k-mer lazy-rc = +120% (2.2×)** (ecoli100 k31, 60 M kmers), **byte-identical ids** (validated vs the per-k-mer path). Now the default `--kmer --kmer-bv` path (no env gate). Note: it emits `pos:presence(1):id[:f|r]`; occurrence counts remain on `--kmer-count --kmer-bv`. Speed summary (canonical MPHF id): 0.48 (always-rc) → 0.75 (lazy-rc) → **1.61 (keep-going)** M/s; ref presence 10, count-bv 2.9, non-canonical-id 1.0. **The fast path is NOT a small tweak:** you can't get `rc(x)` cheaply from a finished fw-only search, so it requires the fw search ITSELF to be bidirectional (`extend_bidirectional` maintains both `fw_interval` and `rc_interval`; then canonical id = `rank(min(lb_fw,lb_rc))`, orientation = `lb_fw<=lb_rc`). That means **restructuring `--kmer-bv` off `query_kmers_from(single=true)` onto a bidirectional walk**. The fast keep-going walks DON'T maintain rc: `query_kmers_count_bv` is fw-only on a *subset* interval, and the existing `query_kmers_from_bidirectional` drops rc on the `partial_matches` left-extension (line ~192 uses `backward_search_step` on fw only; the fix is `extend_bidirectional` there). So the proper speed project = a bidirectional keep-going walk that keeps rc for every k-mer, emitting id+orientation per k-mer. Single-kmer bidirectional (no keep-going) is a smaller intermediate but loses the fw ftab assist. VERIFY rc semantics via the bijection test (already passing on the correctness path, which confirms the `rank(min(lb))` math).
- [x] **Sparse (Elias-Fano) `B_k` for the id** (`build-kmerbv` also writes `kmerbv.<k>.sd`; query selects it via `MOVI_ID_BV=sd`, rank routed through `kmerbv_rank1()`). **Necessary for SSHash size-competitiveness.** ecoli100 k=31 canonical: **sd_vector 32 MB vs dense bv+rank 149 MB (21%), and 0.57× SSHash's 56 MB** (the add-on goes from 2.6× SSHash to *smaller* than SSHash). Ids **identical** to dense (verified); query **~6.5% slower** (sd rank, 1.08→1.01 M/s on the unoptimized single-k-mer path) and **−114 MB RSS** (818→704). canonical mode halves the ones, so sd is even more favorable. (commit pending) — *Note: the on-disk/size-crossover number is now competitive; the remaining ~131 MB `kmerbv_all_p` is RSS-only (run-local-id elimination is a separate follow-up).*
- [ ] `num_kmers` / `print-info` accessor (read from `.meta`).
- [ ] `access`: `select` over `B_k` + forward reconstruction; `access-kmer` CLI.
- [ ] extend the bijection test to canonical mode; cross-check `num_kmers` vs an SSHash `--canonical` index.
- [x] **`kmerbv_all_p` ELIMINATED via stride addressing (BL's idea) — DONE.** Address each `B_k` mark at `run*L + offset` (`L = MAX_RUN_LENGTH = 2047` for regular-thresholds; rows are `<= L` after length-splitting, so the address is injective and order-preserving) instead of its absolute BWT row. The id rank then needs NO per-run position table: `id = rank_sd(run*L+offset + 1) - 1`. Built directly via `sd_vector_builder` over universe `r*L` (monotone addresses); `L` stored in `.meta`; query reads it. **Result (ecoli100 k31): query RSS 704 → 575 MB (−129 MB, all_p gone), speed 1.61 → 1.75 M/s (+9%), ids byte-identical (cross-checked sd==dense).** Tradeoff: the sd universe grows `n → r*L`, so the id structure on disk grows 32 → 62 MB (now ~parity with SSHash 56 MB, was 0.57×); inflation `= L/(n/r)`, so it's small for high-n/r (repetitive) and larger for low-n/r — but where it's larger, `all_p` was even larger, so RSS always drops. Dense rep keeps absolute addressing + all_p (opt-in). The count path was already run-local (unchanged). Possible refinement: use `L = actual max row length` (`<= MAX_RUN_LENGTH`) to cut inflation.
