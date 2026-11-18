@echo off
setlocal EnableExtensions
rem Downloads the memory task recording of every participant: 122 files, about 186 GB.
rem Safe to stop and re-run: finished files are skipped, partial files resume, every file is
rem checked by size and SHA256. At 8 MB/s this takes about 7 hours.
rem Rest and video for everyone (23 GB more):  02b_download_cohort.cmd --tasks rest video
pushd "%~dp0..\.." || exit /b 1
if /i "%HB_CONFIRM_DOWNLOAD%"=="YES" goto :confirmed
echo This step downloads large recordings from OpenNeuro (see the comment at the top of this file).
set /p "ANS=Type YES to start the download: "
if /i not "%ANS%"=="YES" (echo Download not started. & popd & exit /b 4)
:confirmed
call "%~dp0_env.cmd" || goto :fail
"%PY%" -m hb.download_ds008302 --mode cohort --allow-large %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DOWNLOAD COMPLETE AND VERIFIED.
if not "%RC%"=="0" echo NOT COMPLETE, exit code %RC%. Re-run this script to resume.
popd
exit /b %RC%
:fail
popd
exit /b 1
