#!/bin/bash
# install script for lutiffan/methylModes
# usage: install_methylModes.sh [--yes] [--repo user/repo] [--r-min x.y.z]

set -euo pipefail

REPO="lutiffan/methylModes"
R_MIN="4.0.0"
AUTO_YES=0

# Simple output helpers
info(){ printf "\033[1;34m[INFO]\033[0m %s\n" "$*"; }
warn(){ printf "\033[1;33m[WARN]\033[0m %s\n" "$*"; }
err(){ printf "\033[1;31m[ERROR]\033[0m %s\n" "$*"; }
die(){ err "$*"; exit 1; }

usage(){
  cat <<EOF
Usage: $(basename "$0") [--yes] [--repo user/repo] [--r-min x.y.z]

--yes        : auto-install missing R helper packages (remotes) non-interactively
--repo       : GitHub repo (default: ${REPO})
--r-min      : minimum R version required (default: ${R_MIN})
EOF
  exit 1
}

# parse args
while [ $# -gt 0 ]; do
  case "$1" in
    --yes) AUTO_YES=1; shift;;
    --repo) REPO="$2"; shift 2;;
    --r-min) R_MIN="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

# Enforce minimum major-4 baseline: require at least 4.0.0
MIN_MAJOR4="4.0.0"
ver_ge() {
  # compare versions, return 0 if $1 >= $2
  [ "$(printf '%s\n' "$1" "$2" | sort -V | head -n1)" != "$1" ] || [ "$1" = "$2" ]
}
if ! ver_ge "${R_MIN}" "${MIN_MAJOR4}"; then
  warn "Requested minimum R version (${R_MIN}) is below ${MIN_MAJOR4}; elevating required R version to ${MIN_MAJOR4}."
  R_MIN="${MIN_MAJOR4}"
fi

# 1) Check R availability and version
if ! command -v R >/dev/null 2>&1 && ! command -v Rscript >/dev/null 2>&1; then
  die "R not found in PATH. Please install R (>= ${R_MIN})."
fi

# Prefer 'R' binary and parse its --version output; fall back to Rscript --version
if command -v R >/dev/null 2>&1; then
  R_VER_RAW=$(R --version 2>&1 | awk '/^R version/ {print $3; exit}' || true)
else
  # fallback: Rscript --version outputs something like:
  # "R scripting front-end version 4.2.2 (2022-10-31)"
  R_VER_RAW=$(Rscript --version 2>&1 | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?' | head -n1 || true)
fi

if [ -z "${R_VER_RAW}" ]; then
  die "Unable to detect R version. Ensure 'R' or 'Rscript' is in PATH."
fi

info "Detected R version: ${R_VER_RAW}"
if ! ver_ge "${R_VER_RAW}" "${R_MIN}"; then
  die "R >= ${R_MIN} required. Please upgrade R."
fi

# 2) Check compilers / build tools (prefer GNU, allow clang family fallback)
MISSING_BUILD_TOOLS=()
# preferred (GNU) and fallbacks (clang family)
CC_BIN=""
CXX_BIN=""
FC_BIN=""

if command -v gcc >/dev/null 2>&1; then
  CC_BIN="gcc"
elif command -v clang >/dev/null 2>&1; then
  CC_BIN="clang"
fi

if command -v g++ >/dev/null 2>&1; then
  CXX_BIN="g++"
elif command -v clang++ >/dev/null 2>&1; then
  CXX_BIN="clang++"
fi

if command -v gfortran >/dev/null 2>&1; then
  FC_BIN="gfortran"
elif command -v flang >/dev/null 2>&1; then
  FC_BIN="flang"
fi

# collect missing items
[ -z "$CC_BIN" ] && MISSING_BUILD_TOOLS+=("C compiler (gcc/clang)")
[ -z "$CXX_BIN" ] && MISSING_BUILD_TOOLS+=("C++ compiler (g++/clang++)")
[ -z "$FC_BIN" ] && MISSING_BUILD_TOOLS+=("Fortran compiler (gfortran/flang)")

if [ ${#MISSING_BUILD_TOOLS[@]} -ne 0 ]; then
  warn "Missing build tools: ${MISSING_BUILD_TOOLS[*]}"
  warn "These are required to build R packages from source (common on Linux)."
  if command -v apt-get >/dev/null 2>&1; then
    info "On Debian/Ubuntu you can install them with (GNU toolchain recommended):"
    printf "  sudo apt-get update && sudo apt-get install -y build-essential gfortran libcurl4-openssl-dev libxml2-dev libssl-dev\n"
    printf "  (or to use clang: sudo apt-get install -y clang clang++ gfortran ...)\n"
  elif command -v yum >/dev/null 2>&1; then
    info "On CentOS/RHEL you can install them with:"
    printf "  sudo yum groupinstall -y 'Development Tools' && sudo yum install -y gcc-gfortran libcurl-devel libxml2-devel openssl-devel\n"
    printf "  (or to use clang: sudo yum install -y clang clang-devel gfortran ...)\n"
  else
    info "Please install a C/C++/Fortran toolchain and the development headers for curl/xml/ssl."
  fi
fi

# If using clang-family fallbacks, export CC/CXX/FC so R's package build picks them up
if [ -n "$CC_BIN" ] && [ -n "$CXX_BIN" ]; then
  # only export if not using default GNU names (user may still want system defaults)
  if [[ "$CC_BIN" == clang* || "$CXX_BIN" == clang* ]]; then
    export CC="$CC_BIN"
    export CXX="$CXX_BIN"
    info "Using Clang family compilers: CC=$CC CXX=$CXX (exported for R builds)"
  else
    info "Using GNU compilers: CC=$CC_BIN CXX=$CXX_BIN"
  fi
fi

if [ -n "$FC_BIN" ]; then
  export FC="$FC_BIN"
  info "Fortran compiler detected: FC=$FC (exported)"
else
  warn "No Fortran compiler detected. Some R packages require gfortran/flang to build native code."
fi

# 3) Check git / network
if ! command -v git >/dev/null 2>&1; then
  warn "git not found. remotes::install_github can still use the tarball, but git is recommended."
fi

# 4) Check R helper package (remotes preferred; devtools acceptable)
R_CHECK_REPOS_SCRIPT=$(cat <<'RRC'
has_remotes <- requireNamespace("remotes", quietly=TRUE)
has_devtools <- requireNamespace("devtools", quietly=TRUE)
if (has_remotes && has_devtools) cat("BOTH")
else if (has_remotes) cat("REMOTES_ONLY")
else if (has_devtools) cat("DEVTOOLS_ONLY")
else cat("NONE")
RRC
)

R_CHK_OUT=$(Rscript -e "$R_CHECK_REPOS_SCRIPT" 2>/dev/null || true)

case "$R_CHK_OUT" in
  BOTH)
    info "Found both remotes and devtools in R environment."
    ;;
  REMOTES_ONLY)
    info "Found remotes but devtools is missing."
    if [ "$AUTO_YES" -eq 1 ]; then
      if r_install_and_check "devtools"; then
        info "'devtools' installed."
      else
        warn "Automatic install/verify of 'devtools' failed. You can install it in R with: install.packages('devtools')"
      fi
    else
      warn "Run this script with --yes to attempt automatic installation of 'devtools', or in R run:"
      printf "  install.packages('devtools', repos='https://cloud.r-project.org')\n"
    fi
    ;;
  DEVTOOLS_ONLY)
    info "Found devtools but remotes is missing."
    if [ "$AUTO_YES" -eq 1 ]; then
      if r_install_and_check "remotes"; then
        info "'remotes' installed."
      else
        warn "Automatic install/verify of 'remotes' failed. You can install it in R with: install.packages('remotes')"
      fi
    else
      warn "Run this script with --yes to attempt automatic installation of 'remotes', or in R run:"
      printf "  install.packages('remotes', repos='https://cloud.r-project.org')\n"
    fi
    ;;
  NONE|*)
    warn "Neither remotes nor devtools R packages detected."
    if [ "$AUTO_YES" -eq 1 ]; then
      # try remotes first (lighter), then devtools
      ok_any=0
      if r_install_and_check "remotes"; then ok_any=1; fi
      if r_install_and_check "devtools"; then ok_any=1; fi
      if [ "$ok_any" -eq 1 ]; then
        info "At least one installer (remotes/devtools) is now available."
      else
        warn "Automatic installation attempts for remotes/devtools failed. You may need system dependencies for devtools."
      fi
    else
      warn "Run this script with --yes to allow automatic installation of 'remotes' and 'devtools', or in R run:"
      printf "  install.packages(c('remotes','devtools'), repos='https://cloud.r-project.org')\n"
      die "Missing R helper packages 'remotes' and 'devtools'. Aborting."
    fi
    ;;
esac

# Re-check that at least one installer package exists; if not, attempt fallback install from GitHub tarball
R_FINAL_CHECK=$(Rscript -e 'if (requireNamespace("remotes", quietly=TRUE) || requireNamespace("devtools", quietly=TRUE)) cat("OK") else cat("MISSING")' 2>/dev/null || true)
if [ "$R_FINAL_CHECK" != "OK" ]; then
  warn "Neither 'remotes' nor 'devtools' are available in R after attempted installs."
  # If --yes was given, attempt to install remotes first (less heavy than devtools)
  if [ "$AUTO_YES" -eq 1 ]; then
    if r_install_and_check "remotes"; then
      R_FINAL_CHECK="OK"
    else
      warn "Automatic install/verify of 'remotes' failed."
    fi
  fi
fi

# If we still don't have remotes/devtools, fallback to direct tarball install from GitHub
if [ "$R_FINAL_CHECK" != "OK" ]; then
  warn "Falling back to direct GitHub source install (no remotes/devtools)."
  TMPDIR=$(mktemp -d)
  repo_name="${REPO#*/}"
  TAR_URL="https://github.com/${REPO}/archive/refs/heads/main.tar.gz"
  info "Downloading ${TAR_URL} ..."
  if command -v curl >/dev/null 2>&1; then
    curl -L -o "${TMPDIR}/repo.tar.gz" "$TAR_URL" || { warn "curl failed"; }
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${TMPDIR}/repo.tar.gz" "$TAR_URL" || { warn "wget failed"; }
  else
    die "Neither curl nor wget available to download GitHub tarball. Install one or install remotes/devtools in R."
  fi

  if [ ! -s "${TMPDIR}/repo.tar.gz" ]; then
    rm -rf "${TMPDIR}"
    die "Download failed or produced empty file. Cannot proceed with fallback."
  fi

  info "Unpacking archive..."
  tar -xzf "${TMPDIR}/repo.tar.gz" -C "$TMPDIR" || { rm -rf "${TMPDIR}"; die "Failed to unpack archive."; }

  # locate extracted folder (github archive names it as repo-main or repo-<sha>)
  EXDIR=$(find "$TMPDIR" -maxdepth 1 -type d -name "${repo_name}-*" | head -n1 || true)
  if [ -z "$EXDIR" ]; then
    rm -rf "${TMPDIR}"
    die "Could not locate extracted package directory in ${TMPDIR}."
  fi

  info "Running R CMD INSTALL on ${EXDIR} ..."
  if R CMD INSTALL "$EXDIR" >/dev/null 2>&1; then
    info "R CMD INSTALL completed. Verifying by loading the package..."
    if Rscript -e 'if (!requireNamespace("methylModes", quietly=TRUE)) { cat("LOAD_FAIL\n"); quit(status=2) } else { suppressPackageStartupMessages(library(methylModes)); cat("LOADED\n") }'; then
      info "methylModes installed and verified via tarball fallback."
      rm -rf "${TMPDIR}"
      exit 0
    else
      rm -rf "${TMPDIR}"
      die "Installed from tarball but package failed to load. Inspect R output above."
    fi
  else
    rm -rf "${TMPDIR}"
    die "R CMD INSTALL from tarball failed. You may need build dependencies or to install remotes/devtools in R."
  fi
fi

# 5) Install the package from GitHub
info "Installing ${REPO} from GitHub (this may take some time)..."

# Replace heredoc->variable + Rscript -e (can hit quoting/expansion issues) with a temp file.
TMP_R_SCRIPT="$(mktemp /tmp/install_methylModes.R.XXXXXX)"
cat > "$TMP_R_SCRIPT" <<RINS
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    remfun <- devtools::install_github
  } else {
    stop("remotes/devtools missing")
  }
} else {
  remfun <- remotes::install_github
}
# expand ${REPO} here in the heredoc so it's explicit in the script file
remfun("${REPO}", dependencies = TRUE, upgrade = "never")
cat("SUCCESS\n")
RINS

if Rscript "$TMP_R_SCRIPT"; then
  info "Install command completed. Verifying by attempting to load the package in R..."
  # verification: attempt to load package and print version; non-zero exit on failure
  if Rscript -e 'if (!requireNamespace("methylModes", quietly=TRUE)) { cat("LOAD_FAIL: package not found\n"); quit(status=2) } else { suppressPackageStartupMessages(library(methylModes)); cat("LOADED: ", as.character(packageVersion("methylModes")), "\n") }'; then
    info "methylModes installed and verified successfully."
    rm -f "$TMP_R_SCRIPT"
    exit 0
  else
    rm -f "$TMP_R_SCRIPT"
    die "Installation completed but failed to load methylModes in R. Inspect R output above for details."
  fi
else
  rm -f "$TMP_R_SCRIPT"
  die "Installation failed. Inspect the output above to see which dependency or compile step failed."
fi
    exit 0
  else
    rm -f "$TMP_R_SCRIPT"
    die "Installation completed but failed to load methylModes in R. Inspect R output above for details."
  fi
else
  rm -f "$TMP_R_SCRIPT"
  die "Installation failed. Inspect the output above to see which dependency or compile step failed."
fi
