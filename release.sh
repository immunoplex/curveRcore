#!/usr/bin/env bash
# =============================================================================
# release.sh — cut a versioned, tagged release of an R package (curveRcore).
#
# Usage (from the package root, in Git Bash / WSL / macOS / Linux):
#   ./release.sh              # releases version 0.4.0 (default below)
#   ./release.sh 0.4.1        # releases a different version
#   NO_NEWS=1 ./release.sh    # skip the NEWS.md edit
#   NO_PUSH=1 ./release.sh    # do everything locally but don't push / release
#
# What it does, in order:
#   1. Preflight: confirm we're in the curveRcore git repo, tree is clean-ish,
#      and the tag doesn't already exist.
#   2. Bump DESCRIPTION Version to $VERSION.
#   3. Prepend a NEWS.md section (unless NO_NEWS=1).
#   4. Show you the diff and PAUSE for confirmation.
#   5. Commit, create an annotated tag v$VERSION.
#   6. Push the branch and the tag (unless NO_PUSH=1).
#   7. Create a GitHub Release from the NEWS notes if `gh` is installed &
#      authenticated; otherwise print the URL to do it by hand.
#
# Safe to re-read before running. It never force-pushes and never deletes tags.
# =============================================================================
set -euo pipefail

VERSION="${1:-0.4.0}"
TAG="v${VERSION}"
PKG_EXPECTED="curveRcore"

# ── 1. Preflight ─────────────────────────────────────────────────────────────
command -v git >/dev/null || { echo "ERROR: git not found on PATH."; exit 1; }
git rev-parse --is-inside-work-tree >/dev/null 2>&1 \
  || { echo "ERROR: not inside a git repository. cd to the curveRcore repo first."; exit 1; }

# Move to repo root so relative paths (DESCRIPTION, NEWS.md) are correct.
cd "$(git rev-parse --show-toplevel)"

[ -f DESCRIPTION ] || { echo "ERROR: no DESCRIPTION here — is this the package root?"; exit 1; }

PKG_NAME="$(awk -F': *' '/^Package:/{print $2; exit}' DESCRIPTION | tr -d '\r')"
if [ "$PKG_NAME" != "$PKG_EXPECTED" ]; then
  echo "WARNING: DESCRIPTION Package is '$PKG_NAME', expected '$PKG_EXPECTED'."
  read -r -p "Continue anyway? [y/N] " ans; [ "${ans:-N}" = "y" ] || exit 1
fi

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
echo "Repo:    $(git rev-parse --show-toplevel)"
echo "Package: $PKG_NAME"
echo "Branch:  $BRANCH"
echo "Version: $VERSION   Tag: $TAG"
echo

# Tag must not already exist (locally or on the remote).
if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null; then
  echo "ERROR: tag $TAG already exists locally. Aborting."; exit 1
fi
if git ls-remote --exit-code --tags origin "$TAG" >/dev/null 2>&1; then
  echo "ERROR: tag $TAG already exists on origin. Aborting."; exit 1
fi

# ── 2. Bump DESCRIPTION Version ──────────────────────────────────────────────
# Prefer the `desc` package (preserves DCF formatting → clean diff). It ships
# with devtools/usethis, so it's almost certainly installed. Fall back to a
# portable awk rewrite if not.
if command -v Rscript >/dev/null && \
   Rscript -e 'quit(status = as.integer(!requireNamespace("desc", quietly=TRUE)))' >/dev/null 2>&1; then
  Rscript -e "desc::desc_set_version('${VERSION}'); cat('DESCRIPTION Version ->', as.character(desc::desc_get_version()), '\n')"
else
  echo "(desc package not found — editing DESCRIPTION with awk)"
  awk -v v="$VERSION" 'BEGIN{done=0}
    /^Version:/ && !done {print "Version: " v; done=1; next} {print}' \
    DESCRIPTION > DESCRIPTION.tmp && mv DESCRIPTION.tmp DESCRIPTION
fi

# Sanity: DESCRIPTION really says $VERSION now.
DESC_VER="$(awk -F': *' '/^Version:/{print $2; exit}' DESCRIPTION | tr -d '\r')"
[ "$DESC_VER" = "$VERSION" ] || { echo "ERROR: DESCRIPTION Version is '$DESC_VER', expected '$VERSION'."; exit 1; }

# ── 3. NEWS.md ───────────────────────────────────────────────────────────────
NEWS_BODY="$(cat <<EOF
* **Fixed the inflection-point calculation.** \`inflect_x\` no longer snaps to the
  flat lower asymptote (~log10(grid_min_conc) = -4). Added \`compute_inflection()\`,
  which returns the exact closed-form inflection from the fitted parameters
  (x = c for logistic4/gompertz4/loglogistic5; x = c + b*ln(g) for logistic5),
  with y evaluated through the model's own forward function.
* Hardened the grid-based inflection fallback to locate the steepest point
  (argmax |dy/dx|) instead of min |d2y/dx2|, making it immune to flat grid tails.
EOF
)"

if [ "${NO_NEWS:-0}" != "1" ]; then
  DATE="$(date +%Y-%m-%d)"
  HEADER="# ${PKG_NAME} ${VERSION} (${DATE})"
  if [ -f NEWS.md ] && grep -q "^# ${PKG_NAME} ${VERSION}\b" NEWS.md; then
    echo "(NEWS.md already has a ${VERSION} section — leaving it as is)"
  else
    { printf '%s\n\n%s\n\n' "$HEADER" "$NEWS_BODY"; [ -f NEWS.md ] && cat NEWS.md; } > NEWS.md.tmp \
      && mv NEWS.md.tmp NEWS.md
    echo "Prepended ${VERSION} section to NEWS.md"
  fi
fi

# ── 4. Review & confirm ──────────────────────────────────────────────────────
git add -A
echo
echo "──────── staged for release ────────"
git --no-pager status --short
echo "────────────────────────────────────"
git --no-pager diff --cached -- DESCRIPTION NEWS.md || true
echo
read -r -p "Commit, tag ${TAG}, and push? [y/N] " ans
[ "${ans:-N}" = "y" ] || { echo "Aborted before commit. (DESCRIPTION/NEWS edits are staged but uncommitted.)"; exit 1; }

# ── 5. Commit + annotated tag ────────────────────────────────────────────────
git commit -m "Release ${TAG}: analytic inflection point (compute_inflection) + hardened grid fallback"
git tag -a "$TAG" -m "${PKG_NAME} ${VERSION}

${NEWS_BODY}"
echo "Created annotated tag ${TAG}."

# ── 6. Push ──────────────────────────────────────────────────────────────────
if [ "${NO_PUSH:-0}" = "1" ]; then
  echo "NO_PUSH=1 set — skipping push. To push later:"
  echo "  git push origin ${BRANCH} && git push origin ${TAG}"
  exit 0
fi
git push origin "$BRANCH"
git push origin "$TAG"
echo "Pushed ${BRANCH} and ${TAG} to origin."

# ── 7. GitHub Release (optional, via gh CLI) ─────────────────────────────────
if command -v gh >/dev/null && gh auth status >/dev/null 2>&1; then
  printf '%s\n' "$NEWS_BODY" > .release_notes.tmp
  gh release create "$TAG" \
    --title "${PKG_NAME} ${VERSION}" \
    --notes-file .release_notes.tmp \
    --target "$BRANCH"
  rm -f .release_notes.tmp
  echo "GitHub Release ${TAG} created."
else
  REMOTE_URL="$(git remote get-url origin 2>/dev/null || echo '')"
  echo "gh CLI not available/authenticated — the tag is pushed, so create the"
  echo "Release in the browser:"
  case "$REMOTE_URL" in
    git@github.com:*) SLUG="${REMOTE_URL#git@github.com:}"; SLUG="${SLUG%.git}";;
    https://github.com/*) SLUG="${REMOTE_URL#https://github.com/}"; SLUG="${SLUG%.git}";;
    *) SLUG="";;
  esac
  [ -n "$SLUG" ] && echo "  https://github.com/${SLUG}/releases/new?tag=${TAG}"
fi
