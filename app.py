"""Local web server for OligoVisualizer V1: serves the page and the drawing endpoint."""

import socket
import webbrowser

from flask import Flask, jsonify, request

from notation import parse, summary_line
from render import legend_items, render_svg

app = Flask(__name__, static_folder="static", static_url_path="/static")


def _problem_json(problem) -> dict:
    return {
        "message": problem.message,
        "position": problem.position,
        "token": problem.token,
    }


def build_payload(raw: str) -> dict:
    """Turn a raw IDT string into the JSON payload the page consumes."""
    if raw == "":
        return {
            "ok": False,
            "svg": "",
            "summary": "",
            "legend": [],
            "warnings": [],
            "errors": [
                {"message": "no sequence supplied", "position": 0, "token": ""},
            ],
        }

    seq = parse(raw)
    warnings = [_problem_json(problem) for problem in seq.warnings]
    errors = [_problem_json(problem) for problem in seq.errors]

    if not seq.nucleotides:
        if not errors:
            errors = [
                {"message": "no nucleotide found", "position": 0, "token": ""},
            ]
        return {
            "ok": False,
            "svg": "",
            "summary": "",
            "legend": [],
            "warnings": warnings,
            "errors": errors,
        }

    return {
        "ok": True,
        "svg": render_svg(seq),
        "summary": summary_line(seq),
        "legend": legend_items(seq),
        "warnings": warnings,
        "errors": errors,
    }


@app.get("/")
def index():
    """Serve the viewer page verbatim; it is static and never templated."""
    return app.send_static_file("index.html")


@app.get("/api/render")
def api_render():
    """Return the drawing payload for one pasted sequence; always HTTP 200."""
    return jsonify(build_payload(request.args.get("seq", "")))


def main() -> None:
    """Take the first free port in 8765-8785, announce it, open the browser, serve."""
    port = None
    for candidate in range(8765, 8786):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
                probe.bind(("127.0.0.1", candidate))
            port = candidate
            break
        except OSError:
            continue

    if port is None:
        print("No free port found between 8765 and 8785.")
        return

    url = f"http://127.0.0.1:{port}/"
    print(f"Oligo viewer running at {url}")
    print("Close this window to stop the viewer.")
    webbrowser.open(url)
    app.run(host="127.0.0.1", port=port, debug=False, use_reloader=False)


if __name__ == "__main__":
    main()
