@echo off
setlocal EnableExtensions
rem Creates .venv in the release folder and installs the packages.
rem requirements-lock.txt (exact versions) is used only when this Python has the same minor version
rem as the one that built the lock file; otherwise requirements.txt (lower bounds) is used.
rem To choose the interpreter: set HB_PYTHON=C:\path\to\python.exe before running.
rem To skip this step and use an environment you already have: set HB_PYTHON to its python.exe and
rem run scripts\run_tests.cmd and scripts\verify_release.cmd directly.
pushd "%~dp0.." || exit /b 1
set "BASE=python"
if defined HB_PYTHON set "BASE=%HB_PYTHON%"
if not exist ".venv\Scripts\python.exe" (
  "%BASE%" -m venv .venv || (echo ERROR: could not create .venv with "%BASE%". Install Python 3.10 or newer, or set HB_PYTHON. & popd & exit /b 1)
)
".venv\Scripts\python.exe" -m pip install --upgrade pip
set "REQ=requirements.txt"
".venv\Scripts\python.exe" -c "import re,sys;t=open('requirements-lock.txt',encoding='utf-8').read();m=re.search(r'Python (\d+\.\d+)',t);v='{}.{}'.format(*sys.version_info[:2]);ok=bool(m) and m.group(1)==v and re.search(r'(?m)^[A-Za-z0-9_.-]+==',t) is not None;print('Python',v,'; lock file built with Python',m.group(1) if m else 'unknown');sys.exit(0 if ok else 3)"
if "%ERRORLEVEL%"=="0" set "REQ=requirements-lock.txt"
if "%REQ%"=="requirements.txt" echo The lock file does not match this Python; installing lower bounds from requirements.txt instead.
echo Installing from %REQ%
".venv\Scripts\python.exe" -m pip install -r "%REQ%"
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Next: scripts\run_tests.cmd and scripts\verify_release.cmd
if not "%RC%"=="0" echo INSTALL FAILED, exit code %RC%.
popd
exit /b %RC%
