import requests
res = requests.get('https://lasp.colorado.edu/mms/sdc/public/data/sdc/burst/all_mms1_summ/2015/09/01/burst_all_mms1_summ_20150901_121114.png')
type(res)
res.status_code == requests.codes.ok
len(res.text)
print(res.text[:250])
