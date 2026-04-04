#!/bin/bash
set -euo pipefail

# Only run in remote (Claude Code on the web) environments
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# This is a static HTML app with no dependencies to install.
# D3.js is loaded from CDN in index.html.
# Start a local HTTP server so the app can be previewed.
if ! pgrep -f "python3 -m http.server 8080" > /dev/null 2>&1; then
  cd "$CLAUDE_PROJECT_DIR"
  python3 -m http.server 8080 > /tmp/bunnies-server.log 2>&1 &
  echo "Started HTTP server on port 8080"
fi
