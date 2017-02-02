# -*- coding: robot -*-
*** Settings ***
Documentation	blah blah blah
...		blah blah blah
Resource	resource.robot

*** Test Cases ***
Genome ID Search
	[Tags]	public
	Create Session	coge	${API_URL}  verify=true
	${resp}=        Get Request	    coge	 ${SEARCH}/16911
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		results
	Dictionary Should Contain Item	${resp.json()["results"][0]}		id	16911

Undeleted Unrestricted Notebook Search
	[Tags]	public
	Create Session	coge	${API_URL}  verify=true
	${resp}=        Get Request	    coge	 ${SEARCH}/type::notebook deleted::0 arabidopsis restricted::0
	Should Be Equal As Strings	${resp.status_code}	200
	Dictionary Should Contain Key	${resp.json()}		results
