# -*- coding: robot -*-
*** Settings ***
Library         HttpLibrary.HTTP
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Variables ***

*** Test Cases ***
Export GFF to Data Store
        [Tags]  auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${document0}=   Catenate
        ...     {
        ...             "type": "export_gff",
        ...             "parameters": {
        ...                     "dest_type": "irods",
        ...                     "genome_id": 16911,
        ...                     "overwrite": 1
        ...             }
        ...     }
        Should Be Valid JSON    ${document0}
        ${content0}=    Parse Json      ${document0}
        ${headers0}=    Create Dictionary       Content-Type=application/json
        ${resp0}=       Put Request     coge    ${JOBS}/${AUTH_PARAMS}  data=${content0}        headers=${headers0}
        Should Be Equal As Strings      ${resp0.status_code}    201
        Dictionary Should Contain Item	${resp0.json()}	success	True
        Dictionary Should Contain Key   ${resp0.json()}	id
        Set Suite Variable      ${id0}  ${resp0.json()["id"]}
        Log     ${id0}

Genome Add 1
	[Tags]	auth-required
	Create Session  coge    ${API_URL}  verify=true
        ${content1}=	Evaluate        json.load(open('${DATA_PATH}/genome_add_16911.json', 'r'))       json
	${headers1}=	Create Dictionary	Content-Type=application/json
	${resp1}=	Put Request	coge	${GENOMES}/${AUTH_PARAMS}	data=${content1}	headers=${headers1}
	Should Be Equal As Strings	${resp1.status_code}	201
	Dictionary Should Contain Item  ${resp1.json()}  success	True
	Dictionary Should Contain Key	${resp1.json()}	id
	Set Suite Variable      ${id1}	${resp1.json()["id"]}
        Log     ${id1}

Genome Add 2
        [Tags]  auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${content2}=     Evaluate        json.load(open('${DATA_PATH}/genome_add_16911.json', 'r'))       json
        ${headers2}=     Create Dictionary       Content-Type=application/json
        ${resp2}=        Put Request     coge    ${GENOMES}/${AUTH_PARAMS}       data=${content2}	headers=${headers2}
        Should Be Equal As Strings      ${resp2.status_code}     201
        Dictionary Should Contain Item  ${resp2.json()}  success	True
        Dictionary Should Contain Key   ${resp2.json()}  id
        Set Suite Variable      ${id2}	${resp2.json()["id"]}
        Log     ${id2}

Annotation Add 1
        [Tags]  auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${content3}=    Evaluate        json.load(open('${DATA_PATH}/load_annotation.json', 'r'))       json
        ${headers3}=    Create Dictionary       Content-Type=application/json
	Set To Dictionary	${content3["parameters"]}	genome_id=28114
        ${resp3}=       Put Request     coge	${JOBS}/${AUTH_PARAMS}	data=${content3}	headers=${headers3}
        Should Be Equal As Strings	${resp3.status_code}	201
        Dictionary Should Contain Item	${resp3.json()}	success	True
        Dictionary Should Contain Key	${resp3.json()}	id
        Set Suite Variable	${id3}	${resp3.json()["id"]}
	Log	${id3}

Annotation Add 2
        [Tags]  auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${content4}=    Evaluate        json.load(open('${DATA_PATH}/load_annotation.json', 'r'))       json
        ${headers4}=    Create Dictionary       Content-Type=application/json
        Set To Dictionary       ${content4["parameters"]}       genome_id=28115
        ${resp4}=       Put Request     coge    ${JOBS}/${AUTH_PARAMS}  data=${content4}        headers=${headers4}
        Should Be Equal As Strings      ${resp4.status_code}    201
        Dictionary Should Contain Item	${resp4.json()}	success	True
        Dictionary Should Contain Key   ${resp4.json()}	id
        Set Suite Variable      ${id4}  ${resp4.json()["id"]}
        Log     ${id4}
