@echo off
setlocal EnableExtensions
rem Runs one analysis step of the pipeline from this release folder.
rem   scripts\run_analysis.cmd            lists the steps in order
rem   scripts\run_analysis.cmd 08_run_cohort [extra options passed to the step]
rem Raw data location: set HB_DATA_ROOT=D:\data\ds008302, or config\local_paths.json (see config\example_config.json).
rem Outputs go to derivatives\, results\<stage>\ and logs\ here; results\tables\ (the stored results) is never touched.
if "%~1"=="" goto :list
set "STEP=%~1"
shift
set "ARGS="
:collect
if "%~1"=="" goto :run
set "ARGS=%ARGS% %1"
shift
goto :collect
:run
set "F=%~dp0pipeline\%STEP%.cmd"
if not exist "%F%" (
  echo Unknown step "%STEP%". Run scripts\run_analysis.cmd without arguments to list the steps.
  exit /b 2
)
call "%F%" %ARGS%
exit /b %ERRORLEVEL%
:list
echo Steps, in order (see docs\METHODS.md; times are those of the authors' runs with 6 workers):
echo   00_check_environment          environment and data location
echo   00b_run_tests                 unit tests
echo   01_download_metadata          metadata (small)
echo   02_download_pilot             3 pilot participants (about 5 GB)
echo   02b_download_cohort           cohort selection (185.9 GB; asks nothing, so run it only when you mean to)
echo   03_inventory_pilot, 03b_correct_beats, 05_test_wavelet_core, 05b_null_benchmark, 05c_null_benchmark_real,
echo   06_slow_signals, 07_pilot_coherence   pilot and benchmarks
echo   08_run_cohort                 primary analysis, per participant (hours)
echo   09_cohort_stats               group tests (minutes)
echo   10_robust_reprocess, 11_robust_coherence, 12_null_calibration   sensitivity checks
echo   13_joint_calibration, 14_null_validation                         joint pipeline and synthetic validation
echo   15a_realism_decomposition, 15b_shuffle_validation, 15c_shuffle_real   trial shuffle null (about 8 h)
echo   16a_ica_fit (about 5 h), 16b_ica_apply, 16c_ica_analyse (about 6 h)   ICA cleaning
echo   18a_symmetric_synthetic (about 9 h), 18b_symmetric_real (about 8 h)    symmetric tests
echo   17_figures                    figures of the earlier manuscript version
exit /b 0
