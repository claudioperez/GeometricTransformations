input="$1"
shift
python -m veux.motion $input $input --style frame.surface.scale=10 \
       --recover-rotations iter --show origin.axes $@
