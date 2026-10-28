"""Development launcher with optional auto-login from a gitignored .env.

Normal use is still ``Start Oligo Viewer.bat``.  This launcher exists only for
iterative local testing; it reads IDT_USERNAME and IDT_PASSWORD from ``.env``
and signs in before serving.  It never prints or returns those values.
"""

from __future__ import annotations

from pathlib import Path

import app


def read_env(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            value = value[1:-1]
        values[key] = value
    return values


def main() -> None:
    env = read_env(Path(__file__).with_name(".env"))
    username = env.get("IDT_USERNAME", "")
    password = env.get("IDT_PASSWORD", "")

    print("Development viewer: auto-login is optional; normal use is Start Oligo Viewer.bat.")
    if username and password:
        try:
            app._session.sign_in(username, password)
        except Exception as exc:
            print(f"IDT auto-login failed: {exc}")
        else:
            print("IDT auto-login succeeded.")
        finally:
            env.clear()
            username = ""
            password = ""
    else:
        env.clear()
        print("No IDT_USERNAME/IDT_PASSWORD in .env; use the local sign-in form.")
    print("Development viewer running at http://127.0.0.1:8767/")
    app.app.run(host="127.0.0.1", port=8767, debug=False, use_reloader=False)


if __name__ == "__main__":
    main()
