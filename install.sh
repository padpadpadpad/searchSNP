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
say "Run it now: ${install_dir}/${BIN} --version"

case ":${PATH}:" in
  *":${install_dir}:"*) ;;
  *)
    shell_name="$(basename "${SHELL:-sh}")"
    say ""
    say "${install_dir} is not on your PATH, so '${BIN}' on its own will not be found yet."
    case "$shell_name" in
      fish)
        say "Add it permanently with:"
        say "  fish_add_path ${install_dir}"
        ;;
      zsh)
        say "Add it permanently with:"
        say "  echo 'export PATH=\"${install_dir}:\$PATH\"' >> ~/.zshrc"
        say "  source ~/.zshrc"
        ;;
      bash)
        say "Add it permanently with:"
        say "  echo 'export PATH=\"${install_dir}:\$PATH\"' >> ~/.bashrc"
        say "  source ~/.bashrc"
        ;;
      *)
        say "Add this line to your shell startup file:"
        say "  export PATH=\"${install_dir}:\$PATH\""
        ;;
    esac
    say "Then open a new terminal. Until then, use the full path shown above."
    ;;
esac

"${install_dir}/${BIN}" --version
