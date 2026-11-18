@echo off
setlocal EnableExtensions
rem Creates the project .venv if missing, then checks Python, disk space and
rem connectivity to the two data hosts. Downloads a few kilobytes only.
pushd "%~dp0..\.." || goto :nofolder

if exist ".venv\Scripts\python.exe" goto :havevenv
set "BASEPY="
py -3 --version >nul 2>&1 && set "BASEPY=py -3"
if defined BASEPY goto :mkvenv
python --version >nul 2>&1 && set "BASEPY=python"
if defined BASEPY goto :mkvenv
echo ERROR: no Python 3 found. Install Python 3.10 or newer from python.org, then re-run.
popd
exit /b 1

:mkvenv
echo Creating project virtual environment in .venv using: %BASEPY%
%BASEPY% -m venv ".venv"
if errorlevel 1 goto :venvfail

:havevenv
call "%~dp0_env.cmd" || goto :fail
"%PY%" -m hb.check_environment
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%

:venvfail
echo ERROR: could not create .venv. Copy the full output above and send it back.
:fail
popd
exit /b 1
:nofolder
echo ERROR: cannot enter the project folder.
exit /b 1
