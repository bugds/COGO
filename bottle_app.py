from bottle import Bottle, request, template

from fileWork import (
    makeOptionList,
)

from cogCreator import (
    getSequences,
    getIsoforms
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
        species = set([p.species for p in proteins.values()])
        return template('listbox', options=makeOptionList(species))

    else:
        return 'Provide your e-mail and one of the files below'

@app.get("/calculate")
def calculate():
    ref = request.query.get('hidden')

    if not ref:
        return "Empty params"

    data = {
        'ref': ref
    }

    return template('result', **data)

application = app

# for local usage:
# from bottle import run
# run(app, hostname='localhost', port=9999)