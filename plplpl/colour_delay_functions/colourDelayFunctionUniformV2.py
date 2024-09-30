from plplpl.base_functions import BaseSimpleColourFunction

class ColourDelayFunctionUniformV2(BaseSimpleColourFunction):
    def __init__(self):
        super().__init__("colourDelayFunctionUniformV2", 1, {"colour_min":6, "colour_max":30}, 6, 30, {i:(1/25) for i in range(6, 31)})
