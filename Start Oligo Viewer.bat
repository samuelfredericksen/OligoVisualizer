@echo off
setlocal
cd /d "%~dp0"
title Oligo Viewer

if not exist ".venv\Scripts\python.exe" (
    echo Setting up Python for this folder. This happens once and needs internet.
    py -m venv .venv
    if errorlevel 1 goto :fail
)

".venv\Scripts\python.exe" -c "import flask" >nul 2>&1
if errorlevel 1 (
    echo Installing Flask. This happens once and needs internet.
    ".venv\Scripts\python.exe" -m pip install --disable-pip-version-check -r requirements.txt
    if errorlevel 1 goto :fail
)

".venv\Scripts\python.exe" app.py
goto :end

:fail
echo.
echo Setup failed. Take a screenshot of this window and send it to Sam's assistant.
pause

:end
