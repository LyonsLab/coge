# -*- coding: robot -*-
*** Settings ***
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Variables ***
${EXPERIMENT_ID}	1234

*** Test Cases ***
Experiment Search
	Create Session	coge	${API_URL}  verify=true
	${resp}=        Get Request	    coge	 ${EXPERIMENTS}/search/1234
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		experiments

Experiment Fetch
	Create Session  coge    ${API_URL}  verify=true
	${resp}=	Get Request	coge	${EXPERIMENTS}/31
	Should Be Equal As Strings	${resp.status_code}	200
	${expected}=	Evaluate	json.load(open('${DATA_PATH}/experiment_fetch.json', 'r'))	json
	Dictionaries Should Be Equal	${resp.json()}	${expected}

Experiment Add
    Create Session  coge    ${API_URL}  verify=true
        ${content}=	Evaluate        json.load(open('${DATA_PATH}/experiment_add.json', 'r'))       json
	${headers}=	Create Dictionary	Content-Type=application/json
	${resp}=	Put Request	coge	${EXPERIMENTS}/${AUTH_PARAMS}	data=${content}	headers=${headers}
	Should Be Equal As Strings	${resp.status_code}	201
	Dictionary Should Contain Item	${resp.json()}	success	True
	Dictionary Should Contain Key	${resp.json()}	id

