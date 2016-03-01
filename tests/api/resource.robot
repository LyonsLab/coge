*** Settings ***
Documentation	Blah blah blah
Library		RequestsLibrary
Library		Collections

*** Variables ***
${AUTH_USERNAME}	mbomhoff
${AUTH_TOKEN}	10f4578e8327edc4ce4bb49215147e85
${AUTH_PARAMS}	?username=${AUTH_USERNAME}&token=${AUTH_TOKEN}
${DATA_PATH}	./data
${API_URL}	https://geco.iplantcollaborative.org/coge/api/v1
${ORGANISMS}	/organisms
${GENOMES}      /genomes
${FEATURES}	/features
${EXPERIMENTS}	/experiments
${NOTEBOOKS}	/notebooks
${USERS}	/users
${GROUPS}	/groups
${JOBS}		/jobs
${IRODS}	/irods
${JBROWSE_CONFIG}	/jbrowse/config/tracks
${JBROWSE_SEQUENCE}	/jbrowse/sequence
${JBROWSE_FEATURE}	/jbrowse/track/annotation
${JBROWSE_EXPERIMENT}	/jbrowse/experiment
${JBROWSE_NOTEBOOK}	/jbrowse/experiment/notebook

*** Keywords ***

