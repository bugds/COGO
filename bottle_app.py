from bottle import Bottle, request, template

from fileWork import (
    makeOptionList,
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
    text = analysis.file.getvalue().decode()
    return template('listbox', options=makeOptionList(text))

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