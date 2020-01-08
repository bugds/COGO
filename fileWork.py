def makeOptionList(someList):
    options = ''
    template = '''
                <option value="{}">{}</option>'''
    for element in someList:
        options += template.format(element, element)
    return options