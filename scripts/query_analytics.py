import httplib2
import pandas as pd

from apiclient.discovery import build
from oauth2client.service_account import ServiceAccountCredentials
from sys import argv

# NOTES:
# pip install --upgrade google-api-python-client


def initialize_analyticsreporting(service_account_file, scopes, discovery_uri):
    """Initializes an analyticsreporting service object.

    Returns:
      analytics an authorized analyticsreporting service object.
    """

    credentials = ServiceAccountCredentials.from_json_keyfile_name(service_account_file, scopes=scopes)
    http = credentials.authorize(httplib2.Http())

    # Build the service object.
    analytics = build('analytics', 'v4', http=http, discoveryServiceUrl=discovery_uri)

    return analytics


def get_report(analytics, view_id, startdate_one, startdate_two):
    """Get analytics report for two date ranges
    Ranges start on "startdate"s (in 'YYYY-MM-DD' or 'XdaysAgo' format) and end 'today'.

    :param analytics: authorized analyticsreporting service object (use initialize_analyticsreporting()).
    :param startdate_one: start date for first date range (use YYYY-MM-DD or XdaysAgo).
    :param startdate_two: start date for second date range (use YYYY-MM-DD or XdaysAgo).
    :return: analytics report.
    """

    return analytics.reports().batchGet(
        body={
            'reportRequests': [
                {
                    'viewId': view_id,
                    'dateRanges': [{'startDate': startdate_one,
                                    'endDate': 'today'},
                                   {'startDate': startdate_two,
                                    'endDate': 'today'},

                                   ],
                    'dimensions': [{'name': 'ga:country'}],
                    'dimensionFilterClauses': [{'filters': [{'dimensionName': 'ga:country',
                                                             'not': True,
                                                             'operator': 'EXACT',
                                                             'expressions': ['(not set)']
                                                             }]}],
                    'metrics': [{'expression': 'ga:sessions'}]
                }
            ]
        }
    ).execute()


def write_results(responses, output_filename):
    """Write analytic report to CSV

    :param response1: list of responses (i.e. [response] or [response1, response2]
    :param response2:
    :param output_filename: filename (include path) to write CSV
    :return: true on success
    """

    # Object to store parsed data.
    data = {}  # { country: [all, year, month, week], ... }

    # Extract data from reports.
    reports = []
    for r in responses:
        reports.append(r.get('reports', [])[0])
    #report1 = response1.get('reports', [])[0]  # contains data: All, Year
    #report2 = response2.get('reports', [])[0]  # contains data: Month, Week
    #reports = [report1, report2]
    for report in reports:
        rows = report.get('data', {}).get('rows', [])
        for row in rows:
            dimensions = row.get('dimensions', [])
            metrics = row.get('metrics', [])

            # Save data for each country
            country = dimensions[0].encode('utf-8')
            if data.get(country) == None:
                #print(country)
                data[country] = []
                for m in metrics:
                    sessions = m.get('values')[0]
                    data[country].append(sessions)
            else:
                #print("Already found...appending.")
                for m in metrics:
                    sessions = m.get('values')[0]
                    data[country].append(sessions)

    # Convert dictionary to a dataframe.
    df_base = pd.DataFrame(data.items(), columns=["Country", "Sessions"])
    df = pd.concat([df_base["Country"], df_base["Sessions"].apply(pd.Series)], axis=1)
    df.columns = ["Country", "All", "Year", "Month", "Week"]
    df.sort_values(by="Country", inplace=True)
    df.fillna(value=0, inplace=True)

    # Write results to CSV.
    df.to_csv(output_filename, sep=',', index=False)

    return True


def main(service_account, output_file):
    """ Main function

    :param filename: filename (include path)
    :return:
    """
    service_account_file = service_account
    scopes = ['https://www.googleapis.com/auth/analytics.readonly']
    discovery_uri = ('https://analyticsreporting.googleapis.com/$discovery/rest')
    view_id = '7550812'

    analytics = initialize_analyticsreporting(service_account_file, scopes, discovery_uri)
    response1 = get_report(analytics, view_id, '2005-01-01', '365daysAgo')
    response2 = get_report(analytics, view_id, '30daysAgo', '7daysAgo')
    write_results([response1, response2], output_file)


if __name__ == '__main__':
    if len(argv) == 3:
        main(argv[1], argv[2])
    else:
        print("Error: Must specify service account info & output filename.")
        print("usage: query_analytics.py <path/to/service-account.json> <path/to/output_filename>")
