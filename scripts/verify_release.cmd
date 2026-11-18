@echo off
setlocal EnableExtensions
rem Verifies the release: required files, checksums, privacy scan, portable paths, and regeneration of
rem Table 3, figure source data and manuscript from the stored tables (byte for byte comparison).
rem Writes results\validation\verify_report_latest.md. Exit code 0 only if every check passes.
call "%~dp0_release_env.cmd"
pushd "%HB_ROOT%" || exit /b 1
"%PY%" -W ignore -m hb.verify_release --root "%HB_ROOT%" %*
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
