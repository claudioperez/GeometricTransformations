input="$1"
shift
python -m veux.motion2 $input $input --style frame.surface.scale=2 --recover conv --show origin.axes $@
