#!/usr/bin/env bash
# Install the latest searchSNP release for macOS or Linux.
#
#   curl -fsSL https://raw.githubusercontent.com/padpadpadpad/searchSNP/main/install.sh | bash
#
# Environment variables:
#   INSTALL_DIR   where to put the binary (default: ~/.local/bin)
set -euo pipefail

REPO="padpadpadpad/searchSNP"
BIN="searchSNP"

say() { printf '%s\n' "$*"; }
die() { printf 'error: %s\n' "$*" >&2; exit 1; }

os="$(uname -s)"
arch="$(uname -m)"

case "$os" in
  Darwin) target_os="apple-darwin" ;;
  Linux)  target_os="unknown-linux-musl" ;;
  *) die "unsupported operating system: $os (prebuilt binaries are for macOS and Linux)" ;;
esac

case "$arch" in
  x86_64|amd64)  target_arch="x86_64" ;;
  arm64|aarch64) target_arch="aarch64" ;;
  *) die "unsupported architecture: $arch" ;;
esac

target="${target_arch}-${target_os}"
url="https://github.com/${REPO}/releases/latest/download/${BIN}-${target}.tar.gz"
install_dir="${INSTALL_DIR:-$HOME/.local/bin}"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

say "Downloading ${url}"
curl -fsSL "$url" -o "$tmp/pkg.tar.gz" \
  || die "download failed; is there a release yet? See https://github.com/${REPO}/releases"

tar -xzf "$tmp/pkg.tar.gz" -C "$tmp" || die "could not unpack the download"

mkdir -p "$install_dir"
install -m 755 "$tmp/$BIN" "$install_dir/$BIN"

say "Installed ${install_dir}/${BIN}"

case ":${PATH}:" in
  *":${install_dir}:"*) ;;
  *)
    say ""
    say "${install_dir} is not on your PATH. Add it with:"
    say "  echo 'export PATH=\"${install_dir}:\$PATH\"' >> ~/.profile"
    say "(use ~/.zshrc instead of ~/.profile if you use zsh)"
    ;;
esac

"${install_dir}/${BIN}" --version
