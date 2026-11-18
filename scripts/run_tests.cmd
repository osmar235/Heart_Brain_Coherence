@echo off
setlocal EnableExtensions
rem Runs the automated tests (synthetic signals only, no raw data, no internet). About 5 to 10 minutes.
call "%~dp0_release_env.cmd"
pushd "%HB_ROOT%" || exit /b 1
"%PY%" -W ignore -m unittest discover -s code\tests -t code\tests
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo ALL TESTS PASSED.
if not "%RC%"=="0" echo SOME TESTS FAILED, exit code %RC%.
popd
exit /b %RC%
