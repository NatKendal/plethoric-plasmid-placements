from plplpl.base_functions import BaseSimpleColourFunction

class ColourDelayFunctionUniformV2(BaseSimpleColourFunction):
    def __init__(self):
        super().__init__("colourDelayFunctionUniformV2", 1, {"colour_min":6, "colour_max":30}, 6, 30, {i:(1/25) for i in range(6, 31)})

class ColourUniformV2L30U120(BaseSimpleColourFunction):
    def __init__(self):
        super().__init__("colourUniformV2L30U120", 1, {"colour_min":6, "colour_max":24}, 6, 24, {i:(1/19) for i in range(6, 25)})

class ColourUniformV2L30U150(BaseSimpleColourFunction):
    def __init__(self):
        super().__init__("colourUniformV2L30U150", 1, {"colour_min":6, "colour_max":30}, 6, 30, {i:(1/25) for i in range(6, 31)})
