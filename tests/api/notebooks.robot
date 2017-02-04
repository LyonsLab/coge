# -*- coding: robot -*-
*** Settings ***
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Variables ***
${NOTEBOOK_ID}	758

*** Test Cases ***
Notebook Search
	Create Session	coge	${API_URL}  verify=true
	${resp}=        Get Request	    coge	 ${NOTEBOOKS}/search/${NOTEBOOK_ID}
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		notebooks

Notebook Fetch
	Create Session  coge    ${API_URL}  verify=true
	${resp}=	Get Request	coge	${NOTEBOOKS}/${NOTEBOOK_ID}
	Log Dictionary  ${resp.json()}	level=DEBUG
	Should Be Equal As Strings	${resp.status_code}	200
	${expected}=	Evaluate	json.load(open('${DATA_PATH}/notebook_fetch.json', 'r'))	json
	Dictionaries Should Be Equal	${resp.json()}	${expected}

Notebook Add
        Create Session  coge    ${API_URL}	verify=true debug=1
        ${content}=	Evaluate        json.load(open('${DATA_PATH}/notebook_add.json', 'r'))       json
	${headers}=	Create Dictionary	Content-Type=application/json
	${resp}=	Put Request	coge	${NOTEBOOKS}/${AUTH_PARAMS}	data=${content}	headers=${headers}
	Should Be Equal As Strings	${resp.status_code}	201
	Dictionary Should Contain Item	${resp.json()}	success	True
	Dictionary Should Contain Key	${resp.json()}	id
	Set Suite Variable	${id}	${resp.json()["id"]}

Notebook Fetch Added
	Create Session  coge    ${API_URL}	verify=true debug=1
	Set Test Message	${id}
	${resp}=        Get Request     coge    ${NOTEBOOKS}/${id}/${AUTH_PARAMS}
	Should Be Equal As Strings      ${resp.status_code}     200
	${content}=     Evaluate        json.load(open('${DATA_PATH}/notebook_add.json', 'r'))       json
	Dictionary Should Contain Sub Dictionary	${resp.json()}	${content["metadata"]}

Notebook Search Added
        Create Session  coge    ${API_URL}  verify=true
        ${resp}=        Get Request         coge         ${NOTEBOOKS}/search/${id}/${AUTH_PARAMS}
        Should Be Equal As Strings      ${resp.status_code}     200
        Dictionary Should Contain Key   ${resp.json()}          notebooks
