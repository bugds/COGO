from bottle import Bottle, request, template

from fileWork import (
    pickleFile,
    unpickleFile,
    makeOptionList
)

from cogCreator import (
    getSequences,
    getIsoforms,
    countGenes,
    blastAndCalc
)

app = Bottle(__name__)

@app.get("/")
def index():
    return template('index')

@app.post("/reference")
def reference():
    eMail = request.POST.get('eMail')
    analysis = request.files.analysis
    reanalysis = request.files.reanalysis

    if reanalysis and eMail:
        pass
        return 'Work in progress'

    elif analysis and eMail:
        proteins = dict()
        text = analysis.file.getvalue().decode()

        proteins = getSequences(text, proteins)
        proteins = getIsoforms(proteins)
        genesDict = countGenes(proteins)


        data = {
            'eMail': eMail,
            'proteins': pickleFile(proteins),
            'options': makeOptionList(genesDict)
        }

        return template('listbox', **data)

    else:
        return 'Provide your e-mail and one of the files below'

@app.post("/calculate")
def calculate():
    referencial = (request.forms.get('referencial')).split(';')
    eMail = request.forms.get('eMail')
    proteins = unpickleFile(request.forms.get('proteins'))
    blastDict = None

    result = blastAndCalc(referencial, eMail, proteins, blastDict)

    data = {
        'result': result
    }

    return template('result', **data)

application = app

# for local usage:
# from bottle import run
# run(app, hostname='localhost', port=9999)