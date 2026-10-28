# Draw the molecule on the server, not in the browser

The viewer's job is a schematic, and iteration 2 will bolt thermodynamic predictions onto the
same picture. Rendering the SVG in Python keeps parsing, geometry and (later) prediction in one
language with one test suite, and makes the drawing assertable with `unittest` instead of
screenshots. The cost is a round trip per keystroke and a server process that must be running —
accepted, because the server is already required for comparison and the round trip is on
localhost.
