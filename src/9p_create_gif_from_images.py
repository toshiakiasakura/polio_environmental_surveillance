# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import glob
from PIL import Image

# +
fnames= [str(i).rjust(4, "0") + ".png" for i in range(1, 21)]

frames = [Image.open("../dt_tmp/gif_map/" + f) for f in fnames ]
frame_one = frames[0]
frame_one.save("../dt_tmp/gif_map/gif_map_pc25.gif", format="GIF", append_images=frames, 
               save_all=True, duration=400, loop=10
              )

# +
fnames= [str(i).rjust(4, "0") + ".png" for i in range(1, 16)]

frames = [Image.open("../dt_tmp/gif_map/" + f) for f in fnames ]
frame_one = frames[0]
frame_one.save("../dt_tmp/gif_map/gif_map_pc30.gif", format="GIF", append_images=frames, 
               save_all=True, duration=400, loop=10
              )

# +
fnames= [str(i).rjust(4, "0") + ".png" for i in range(1, 8)]

frames = [Image.open("../dt_tmp/gif_map/" + f) for f in fnames ]
frame_one = frames[0]
frame_one.save("../dt_tmp/gif_map/gif_map_pc50.gif", format="GIF", append_images=frames, 
               save_all=True, duration=400, loop=10
              )

# +
fnames= [str(i).rjust(4, "0") + ".png" for i in range(1, 85)]

frames = [Image.open("../dt_tmp/gif_map/" + f) for f in fnames ]
frame_one = frames[0]
frame_one.save("../dt_tmp/gif_map/gif_map_pc100.gif", format="GIF", append_images=frames, 
               save_all=True, duration=400, loop=10
              )
# -


