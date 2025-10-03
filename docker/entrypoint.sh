#!/bin/bash

# Development container entrypoint
# Provides convenience functions for Movi development

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to build Movi
build_movi() {
    local build_type=${1:-Release}
    local build_dir=${2:-build}
    
    print_info "Building Movi (${build_type}) in ${build_dir}/"
    
    mkdir -p "${build_dir}"
    cd "${build_dir}"
    
    cmake .. -DCMAKE_BUILD_TYPE="${build_type}"
    make -j$(nproc)
    
    print_success "Build completed successfully!"
}

# Function to clean build artifacts
clean_build() {
    local build_dir=${1:-build}
    
    print_info "Cleaning build directory: ${build_dir}/"
    
    if [ -d "${build_dir}" ]; then
        rm -rf "${build_dir}"
        print_success "Build directory cleaned!"
    else
        print_warning "Build directory ${build_dir} does not exist."
    fi
}

# Function to show help
show_help() {
    echo "Movi Development Container"
    echo ""
    echo "Available commands:"
    echo "  build [type] [dir]  - Build Movi (Release|Debug, default: Release)"
    echo "  clean [dir]         - Clean build directory"
    echo "  help                - Show this help"
    echo ""
    echo "Examples:"
    echo "  build Debug build-debug"
    echo "  clean build"
}

# Main entrypoint logic
case "${1:-help}" in
    "build")
        build_movi "$2" "$3"
        ;;
    "clean")
        clean_build "$2"
        ;;
    "help"|"--help"|"-h")
        show_help
        ;;
    *)
        if [ $# -eq 0 ]; then
            # No arguments - start interactive shell
            print_info "Starting Movi development environment..."
            print_info "Use 'help' to see available commands"
            exec /bin/bash
        else
            # Pass through any other commands
            exec "$@"
        fi
        ;;
esac

