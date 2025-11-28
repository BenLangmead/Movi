# Movi Project Makefile
# Provides convenient targets for building, testing, and development

.PHONY: help build build-debug build-release build-release-debug build-profile clean dev dev-build

# Default target
help:
	@echo "Movi Project Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  build          - Build Movi (Release mode)"
	@echo "  build-debug    - Build Movi (Debug mode)"
	@echo "  build-release  - Build Movi (Release mode)"
	@echo "  build-profile  - Build Movi (Profile mode with gprof support)"
	@echo "  clean          - Clean build directories"
	@echo ""
	@echo "Development targets:"
	@echo "  dev            - Start development container"
	@echo "  dev-build      - Build in development container"
	@echo ""
	@echo "Examples:"`
	@echo "  make build-debug"
	@echo "  make dev"

# Local build targets
build: build-release

build-release:
	@echo "Building Movi (Release mode)..."
	@mkdir -p build-release
	@cd build-release && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$$(nproc)
	@echo "Build completed: build-release/"

build-release-debug:
	@echo "Building Movi (ReleaseDebug mode)..."
	@mkdir -p build-release-debug
	@cd build-release-debug && cmake .. -DCMAKE_BUILD_TYPE=ReleaseDebug && make -j$$(nproc)
	@echo "Build completed: build-release-debug/"

build-debug:
	@echo "Building Movi (Debug mode)..."
	@mkdir -p build-debug
	@cd build-debug && cmake .. -DCMAKE_BUILD_TYPE=Debug && make -j$$(nproc)
	@echo "Build completed: build-debug/"

build-profile:
	@echo "Building Movi (Profile mode with gprof support)..."
	@mkdir -p build-profile
	@cd build-profile && cmake .. -DCMAKE_BUILD_TYPE=Profile && make -j$$(nproc)
	@echo "Build completed: build-profile/"

clean:
	@echo "Cleaning build directories..."
	@rm -rf build build-release build-debug build-profile build-release-debug
	@echo "Clean completed"


# Development container targets
dev:
	@echo "Starting development container..."
	@docker build -f docker/Dockerfile.dev .
	@docker run -it --rm \
		-v "$(PWD):/workspace" \
		-w /workspace \
		$(shell docker build -f docker/Dockerfile.dev -q .)

dev-build:
	@echo "Building in development container..."
	@docker build -f docker/Dockerfile.dev .
	@docker run --rm \
		-v "$(PWD):/workspace" \
		-w /workspace \
		--entrypoint /workspace/docker/entrypoint.sh \
		$(shell docker build -f docker/Dockerfile.dev -q .) \
		build

# Utility targets
docker-clean:
	@echo "Cleaning Docker images and containers..."
	@docker system prune -f
	@docker image prune -f

# Show build information
info:
	@echo "Movi Project Information"
	@echo "======================"
	@echo "Source directory: $(PWD)"
	@echo "Available build directories:"
	@ls -la | grep "^d.*build" || echo "  No build directories found"
	@echo ""
	@echo "Docker images:"
	@docker images | grep movi || echo "  No Movi Docker images found"
	@echo ""
	@echo "Docker containers:"
	@docker ps -a | grep movi || echo "  No Movi Docker containers found"

