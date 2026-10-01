#!/usr/bin/env bash
# exit on failure of any command
set -e

echo "=================================== RUN TESTS ==================================="

pushd tests
echo `pwd`
tree -L 2
pytest --capture=fd -v -rxs
popd

echo "=================================== DONE ==================================="

# exit with the exitcode thrown by pytest
exit $?