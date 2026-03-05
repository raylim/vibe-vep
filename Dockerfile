# Build stage: compile vibe-vep with CGO enabled (required for DuckDB + SQLite).
FROM golang:1.24-bookworm AS builder

WORKDIR /build

# Cache module downloads.
COPY go.mod go.sum ./
RUN go mod download

# Build the binary.
COPY . .
RUN CGO_ENABLED=1 go build \
    -ldflags="-s -w" \
    -o /vibe-vep \
    ./cmd/vibe-vep

# Runtime stage: minimal image with C runtime for CGO dependencies.
FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /vibe-vep /usr/local/bin/vibe-vep

# Default data directory.
VOLUME /data
WORKDIR /data

ENTRYPOINT ["vibe-vep"]
