@echo off
rem Shared setup of the pipeline launchers inside the release folder. Not meant to be run alone.
set "HB_ROOT=%~dp0..\.."
for %%I in ("%HB_ROOT%") do set "HB_ROOT=%%~fI"
set "PY=%HB_ROOT%\.venv\Scripts\python.exe"
if not exist "%PY%" if defined HB_PYTHON set "PY=%HB_PYTHON%"
if not exist "%PY%" set "PY=python"
set "PYTHONPATH=%HB_ROOT%\code\src"
set "PYTHONUTF8=1"
set "HB_PROJECT_ROOT=%HB_ROOT%"
if not exist "%HB_ROOT%\logs" mkdir "%HB_ROOT%\logs"
exit /b 0
