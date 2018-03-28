from bokeh.core.properties import List, Instance
from bokeh.models.annotations import Legend, LegendItem
from bokeh.util.compiler import TypeScript

class ToggleLegend(Legend):
    all_items = List(Instance(LegendItem))

    __implementation__ = TypeScript('''\
import {Legend} from "models/annotations/legend"
import * as p from "core/properties"

export class ToggleLegend extends Legend {
    type: "ToggleLegend"
}
ToggleLegend.define({
    all_items: [p.Array, []],
})
''')
