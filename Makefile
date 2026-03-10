.PHONY: build test lint clean install download-testdata docs docs-build wasm wasm-exec parquet-export cudaops build-cuda bench-cuda

# Binary name
BINARY=vibe-vep

# Build flags
LDFLAGS=-ldflags "-s -w"

build:
	go build $(LDFLAGS) -o $(BINARY) ./cmd/vibe-vep

test:
	go test -v -race -cover ./...

lint:
	golangci-lint run ./...

clean:
	rm -f $(BINARY)
	go clean -testcache

install:
	go install $(LDFLAGS) ./cmd/vibe-vep

# Run tests with coverage report
coverage:
	go test -coverprofile=coverage.out ./...
	go tool cover -html=coverage.out -o coverage.html

# Format code
fmt:
	gofmt -s -w .
	goimports -w .

# Download TCGA test data for validation (~1.6GB)
download-testdata:
	./scripts/download_tcga.sh

# Download dependencies
deps:
	go mod download
	go mod tidy

# Export annotations to Parquet (override INPUT and OUTPUT as needed)
INPUT ?= testdata/tcga/chol_tcga_gdc_data_mutations.txt
OUTPUT ?= annotations.parquet
parquet-export: build
	./$(BINARY) export parquet --canonical --pick -o $(OUTPUT) $(INPUT)

# ---- GPU / CUDA targets ----

CUDA_ARCH ?= sm_80
CUDA_HOME ?= /usr/local/cuda

# Compile the CUDA shared library (requires nvcc).
cudaops:
	nvcc -O3 -arch=$(CUDA_ARCH) --compiler-options -fPIC \
		-I$(CUDA_HOME)/include \
		-shared -o internal/cudaops/libcudaops.so \
		internal/cudaops/codon_cuda.cu \
		-L$(CUDA_HOME)/lib64 -lcudart

# Build the main binary with CUDA support (requires cudaops target first).
build-cuda: cudaops
	CGO_ENABLED=1 \
		CGO_CFLAGS="-I$(CUDA_HOME)/include" \
		CGO_LDFLAGS="-L$(CURDIR)/internal/cudaops -lcudaops -Wl,-rpath,$(CURDIR)/internal/cudaops -L$(CUDA_HOME)/lib64 -lcudart" \
		go build -tags cuda $(LDFLAGS) -o $(BINARY)-cuda ./cmd/vibe-vep

# Run the benchmark suite with and without CUDA.
bench-cuda: cudaops
	@echo "=== CPU fallback benchmarks ==="
	GOTOOLCHAIN=go1.24.0 go test ./internal/cudaops/ -bench=. -benchmem -benchtime=3s
	@echo ""
	@echo "=== CUDA GPU benchmarks ==="
	CGO_ENABLED=1 \
		CGO_CFLAGS="-I$(CUDA_HOME)/include" \
		CGO_LDFLAGS="-L$(CURDIR)/internal/cudaops -lcudaops -Wl,-rpath,$(CURDIR)/internal/cudaops -L$(CUDA_HOME)/lib64 -lcudart" \
		LD_LIBRARY_PATH=$(CURDIR)/internal/cudaops:$(CUDA_HOME)/lib64:$(LD_LIBRARY_PATH) \
		GOTOOLCHAIN=go1.24.0 go test -tags cuda ./internal/cudaops/ -bench=. -benchmem -benchtime=3s

# Build WASM binary for interactive tutorial
wasm:
	GOOS=js GOARCH=wasm go build $(LDFLAGS) -o docs/static/wasm/vibe-vep.wasm ./cmd/vibe-vep-wasm

# Copy wasm_exec.js from Go installation
wasm-exec:
	cp "$$(go env GOROOT)/lib/wasm/wasm_exec.js" docs/static/js/wasm_exec.js

# Local docs preview (with drafts)
docs:
	cd docs && hugo server -D

# Build docs for production (includes WASM)
docs-build: wasm
	cd docs && hugo --minify
