@echo off
setlocal EnableExtensions
rem Regenerates Table 3, the claim checks, Figures 1 to 3 (with source data) and the manuscript from the
rem stored result tables in results\tables. Writes to regenerated\ and never overwrites stored files.
call "%~dp0_release_env.cmd"
pushd "%HB_ROOT%" || exit /b 1
set "HB_RESULTS_ROOT=%HB_ROOT%\results\tables"
set "HB_LOGS_ROOT=%HB_ROOT%\results\validation\logs"
set "HB_NO_EXECUTION_LOG=1"
set "HB_OUT_DIR=%HB_ROOT%\regenerated\tables"
"%PY%" -W ignore -m hb.make_tables_three_verdicts || goto :fail
set "HB_OUT_DIR=%HB_ROOT%\regenerated\figures"
"%PY%" -W ignore -m hb.make_figures_three_verdicts || goto :fail
"%PY%" -W ignore -m hb.build_manuscript --src manuscript\manuscript_three_verdicts.src.md --refs manuscript\references.tsv --table3 regenerated\tables\table_3.md --out regenerated\manuscript_three_verdicts.md || goto :fail
echo.
echo DONE. Regenerated files are in regenerated\. scripts\verify_release.cmd compares them with the stored ones.
popd
exit /b 0
:fail
echo FAILED.
popd
exit /b 1
