@echo off
rem Shared setup of the release launchers. Not meant to be run alone.
rem Python: the interpreter named in HB_PYTHON when set, else the .venv of this folder, else python on PATH.
set "HB_ROOT=%~dp0.."
for %%I in ("%HB_ROOT%") do set "HB_ROOT=%%~fI"
set "PY=python"
if exist "%HB_ROOT%\.venv\Scripts\python.exe" set "PY=%HB_ROOT%\.venv\Scripts\python.exe"
if defined HB_PYTHON set "PY=%HB_PYTHON%"
set "PYTHONPATH=%HB_ROOT%\code\src"
set "PYTHONUTF8=1"
set "HB_PROJECT_ROOT=%HB_ROOT%"
echo Python used: %PY%
exit /b 0
