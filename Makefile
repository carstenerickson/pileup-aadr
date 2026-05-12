.PHONY: refresh-chain-sha test lint typecheck install-dev

refresh-chain-sha:
	@if [ ! -f pileup_aadr/data/hg19ToHg38.over.chain.gz ]; then \
		echo "ERROR: pileup_aadr/data/hg19ToHg38.over.chain.gz not found."; \
		echo "Download from UCSC first:"; \
		echo "  curl -L -o pileup_aadr/data/hg19ToHg38.over.chain.gz \\"; \
		echo "    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"; \
		exit 1; \
	fi
	@shasum -a 256 pileup_aadr/data/hg19ToHg38.over.chain.gz \
		| awk '{print $$1}' > pileup_aadr/data/hg19ToHg38.over.chain.gz.sha256
	@echo "Refreshed SHA: $$(cat pileup_aadr/data/hg19ToHg38.over.chain.gz.sha256)"

install-dev:
	pip install -e ".[dev]"

test:
	pytest -v

lint:
	ruff check pileup_aadr tests
	ruff format --check pileup_aadr tests

typecheck:
	mypy pileup_aadr
