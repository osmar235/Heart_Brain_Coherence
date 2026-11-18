@echo off
setlocal EnableExtensions
rem Step 16a: EEG ICA fit stage (extended infomax, Picard) and ICLabel for all participants. Nothing is removed yet.
rem See docs\ANALYSIS_FREEZE_v2_EEG.md. Installs the pinned packages first. Resumable. Optional: --workers N (default 3, about 5 GB each)
pushd "%~dp0..\.." || exit /b 1
call "%~dp0_env.cmd" || goto :fail
echo Installing or checking mne 1.13.2, mne-icalabel 0.10.0, python-picard 0.8.2, onnxruntime 1.31.1 ...
"%PY%" -m pip install --disable-pip-version-check -q mne==1.13.2 mne-icalabel==0.10.0 python-picard==0.8.2 onnxruntime==1.31.1 scikit-learn matplotlib
if errorlevel 1 (
  echo PACKAGE INSTALL FAILED. Send back this window text.
  goto :fail
)
"%PY%" -W ignore -m hb.ica_v2 fit %*
set "RC=%ERRORLEVEL%"
echo.
if "%RC%"=="0" echo DONE. Send back the newest log in the logs folder.
if not "%RC%"=="0" echo FINISHED WITH FAILURES, exit code %RC%. Send back the newest log in the logs folder.
popd
exit /b %RC%
:fail
popd
exit /b 1
