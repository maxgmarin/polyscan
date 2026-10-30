# polyscan — Reference

## Build

Building needs the [Rust toolchain](https://www.rust-lang.org/tools/install) and a C compiler, because some dependencies compile C code (bzip2, xz and zstd for compressed input, libdeflate for bigBed). On macOS, the Xcode command line tools provide one (`xcode-select --install`).

```bash
cargo build --release              # binary at target/release/polyscan
cargo install --locked --path .    # or install to ~/.cargo/bin (usually on PATH)
```

`--locked` builds with the exact dependency versions in `Cargo.lock`, the same ones the tests and releases use.

To install a tagged release without cloning:

```bash
cargo install --locked --git https://github.com/maxgmarin/polyscan --tag v0.4.0
```

To update an install from a clone, run `git pull` and repeat `cargo install --locked --path .`. To remove it, run `cargo uninstall polyscan`.

### Static Linux build

A static Linux binary (no runtime dependencies, e.g. for an HPC cluster) is cross-compiled from macOS with the musl toolchain:
- `brew install FiloSottile/musl-cross/musl-cross`
- `rustup target add x86_64-unknown-linux-musl`

The linker is set in `.cargo/config.toml`. The `CC`/`AR` variables are needed for the C code above:

```bash
CC_x86_64_unknown_linux_musl=x86_64-linux-musl-gcc \
AR_x86_64_unknown_linux_musl=x86_64-linux-musl-ar \
cargo build --release --target x86_64-unknown-linux-musl
```

The binary is at `target/x86_64-unknown-linux-musl/release/polyscan`, and `file` reports it as `static-pie linked`.
