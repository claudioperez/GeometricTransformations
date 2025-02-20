input="$1"
shift
python -m veux.motion $input $input --style frame.surface.scale=10 --recover iter --show origin.axes $@
