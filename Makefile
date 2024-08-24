.PHONY: docker clean

ORG := "logsdonlab"
PROJECT_NAME := "cdr-finder"
VERSION := "latest"
REPO := "${ORG}/${PROJECT_NAME}"
TAG_NAME := "${REPO}:${VERSION}"

clean:
	rm -rf .snakemake/ results/ logs/ benchmarks/

docker:
	snakemake -p -c 1 --sdm apptainer --configfile test/config/config.yaml --containerize | \
	sudo docker build -t "${TAG_NAME}" -f - .

singularity:
	$(MAKE) docker
	sudo apptainer build "${PROJECT_NAME}.sif" "docker-daemon://${TAG_NAME}"

docker-upload:
	sudo docker image push "${TAG_NAME}"
