@echo off
setlocal EnableExtensions
rem Offline unit tests of the downloader. Uses a local test server, no internet.
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
"%PY%" -m unittest discover -s code\tests -v
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
:fail
popd
exit /b 1
