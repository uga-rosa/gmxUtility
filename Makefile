MAKEFILE_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))

.PHONY: install
install:
	ln -snfv ${MAKEFILE_DIR}/src/adsorption.ts ~/.local/bin/adsorption.ts
