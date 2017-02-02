# -*- coding: robot -*-
*** Settings ***
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Variables ***
${GENOME_ID}	16911

*** Test Cases ***
Genome Search
	[Tags]	public
	Create Session	coge	${API_URL}  verify=true
	${resp}=        Get Request	    coge	 ${GENOMES}/search/${GENOME_ID}
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		genomes

Genome Fetch
	[Tags]  public
	Create Session  coge    ${API_URL}  verify=true
	${resp}=	Get Request	coge	${GENOMES}/${GENOME_ID}
	Should Be Equal As Strings	${resp.status_code}	200
	${expected}=	Evaluate	json.load(open('${DATA_PATH}/genome_fetch.json', 'r'))	json
	Dictionary Should Contain Sub Dictionary	${resp.json()}	${expected}

#Genome Fetch Sequence
#        Create Session  coge    ${API_URL} verify=true
#        ${resp}=        Get Request     coge    ${GENOMES}/16911/sequence
#        Should Be Equal As Strings      ${resp.status_code}     200
#       Set Test Message        ${resp.json()}
#        ${expected}=    Evaluate        json.load(open('genome_fetch_sequence.txt', 'r'))
#        Dictionaries Should Be Equal    ${resp}  ${expected}

Genome Add
	[Tags]	auth-required
        Create Session  coge    ${API_URL}  verify=true
        ${content}=	Evaluate        json.load(open('${DATA_PATH}/genome_add.json', 'r'))       json
	${headers}=	Create Dictionary	Content-Type=application/json
	${resp}=	Put Request	coge	${GENOMES}/${AUTH_PARAMS}	data=${content}	headers=${headers}
	Should Be Equal As Strings	${resp.status_code}	201
	Dictionary Should Contain Item  ${resp.json()}  success	True
	Dictionary Should Contain Key	${resp.json()}	id

