from pathlib import Path

import bokeh.models

def load_js_template(file, key):
    fn = Path(file).parent / f'{key}.js'
    return fn.read_text()

def build_js_callback(file, key, format_kwargs=None, args=None):
    if args is None:
        args = {}
    if format_kwargs is None:
        format_kwargs = {}

    code = load_js_template(file, key)
    for format_kwarg, value in format_kwargs.items():
        to_replace = '{' + format_kwarg + '}'
        code = code.replace(to_replace, str(value))

    name = f'{key}_callback'
    callback = bokeh.models.CustomJS(code=code, args=args)
    callback.name = name

    return callback