define install_script
	@cp src/$(1).py $(PREFIX)/$(1)
	@chmod +x $(PREFIX)/$(1)
	@echo "Installed script $(1) in $(PREFIX)."
endef

define uninstall_script
	@rm $(PREFIX)/$(1)
	@echo "Uninstalled script $(1) in $(PREFIX)."
endef

install:
	$(call install_script,sged-concatenate-alignments)
	$(call install_script,sged-create-sequence-index)
	$(call install_script,sged-create-structure-index)
	$(call install_script,sged-disembl2sged)
	$(call install_script,sged-get-all-pairs)
	$(call install_script,sged-group-test-inclusion)
	$(call install_script,sged-group)
	$(call install_script,sged-liftover-index)
	$(call install_script,sged-merge-indexes)
	$(call install_script,sged-merge)
	$(call install_script,sged-paml2sged)
	$(call install_script,sged-randomize-groups)
	$(call install_script,sged-raser2sged)
	$(call install_script,sged-structure-infos)
	$(call install_script,sged-structure-list)
	$(call install_script,sged-summary)
	$(call install_script,sged-translate-coords)
	$(call install_script,sged-ungroup)

uninstall:
	$(call uninstall_script,sged-concatenate-alignments)
	$(call uninstall_script,sged-create-sequence-index)
	$(call uninstall_script,sged-create-structure-index)
	$(call uninstall_script,sged-disembl2sged)
	$(call uninstall_script,sged-get-all-pairs)
	$(call uninstall_script,sged-group-test-inclusion)
	$(call uninstall_script,sged-group)
	$(call uninstall_script,sged-liftover-index)
	$(call uninstall_script,sged-merge-indexes)
	$(call uninstall_script,sged-merge)
	$(call uninstall_script,sged-paml2sged)
	$(call uninstall_script,sged-randomize-groups)
	$(call uninstall_script,sged-raser2sged)
	$(call uninstall_script,sged-structure-infos)
	$(call uninstall_script,sged-structure-list)
	$(call uninstall_script,sged-summary)
	$(call uninstall_script,sged-translate-coords)
	$(call uninstall_script,sged-ungroup)

