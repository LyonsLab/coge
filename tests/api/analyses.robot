# -*- coding: robot -*-
*** Settings ***
Library		HttpLibrary.HTTP
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Variables ***

*** Test Cases ***
Export FASTA to Data Store
	[Tags]  auth-required
	Create Session  coge    ${API_URL}  verify=true
	${document1}=	Catenate
	...	{ 
	...		"type": "export_fasta", 
	...		"parameters": {
	...			"dest_type": "irods", 
	...			"genome_id": 16911, 
	...			"overwrite": 1 
	...		}
	...	}
	Should Be Valid JSON	${document1}
	${content1}=	Parse Json	${document1}	
        ${headers1}=    Create Dictionary       Content-Type=application/json
        ${resp1}=       Put Request     coge	${JOBS}/${AUTH_PARAMS}	data=${content1}	headers=${headers1}
        Should Be Equal As Strings	${resp1.status_code}	201
        Dictionary Should Contain Item	${resp1.json()}	success	True
        Dictionary Should Contain Key	${resp1.json()}	id
        Set Suite Variable	${id1}	${resp1.json()["id"]}
	Log	${id1}

Export GFF to Data Store
        [Tags]  auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${document2}=   Catenate
        ...     {
        ...             "type": "export_gff",
        ...             "parameters": {
        ...                     "dest_type": "irods",
        ...                     "genome_id": 16911,
        ...                     "overwrite": 1
        ...             }
        ...     }
        Should Be Valid JSON    ${document2}
        ${content2}=    Parse Json      ${document2}
        ${headers2}=    Create Dictionary       Content-Type=application/json
        ${resp2}=       Put Request     coge    ${JOBS}/${AUTH_PARAMS}  data=${content2}        headers=${headers2}
        Should Be Equal As Strings      ${resp2.status_code}    201
        Dictionary Should Contain Item	${resp2.json()}	success	True
        Dictionary Should Contain Key   ${resp2.json()}	id
        Set Suite Variable      ${id2}  ${resp2.json()["id"]}
        Log     ${id2}

Export Experiment to Data Store
        [Tags]  auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${document3}=   Catenate
        ...     {
        ...             "type": "export_experiment",
        ...             "parameters": {
        ...                     "dest_type": "irods",
        ...                     "experiment_id": 1234,
        ...                     "overwrite": 1
        ...             }
        ...     }
        Should Be Valid JSON    ${document3}
        ${content3}=    Parse Json      ${document3}
        ${headers3}=    Create Dictionary       Content-Type=application/json
        ${resp3}=       Put Request     coge    ${JOBS}/${AUTH_PARAMS}  data=${content3}        headers=${headers3}
        Should Be Equal As Strings      ${resp3.status_code}    201
        Dictionary Should Contain Item  ${resp3.json()}	success	True
        Dictionary Should Contain Key   ${resp3.json()}	id
        Set Suite Variable      ${id3}  ${resp3.json()["id"]}
        Log     ${id3}
