
% import os

<html>

	<head>
		<script src='static/jquery.js'></script>
		<script src='www/xyz.all.js'></script>
		<script>
		</script>

	</head>

	<body>

		<form>
			Enter element names, atomic symbols, or atomic numbers:
			<input type='text' name='filter' value="{{' '.join([f for f in filter])}}"><br>
		</form>

		<h3>
			{{", ".join(filter)}} 
		</h3>

		<p>
			% for r in results:
				<canvas class='xyz' url='{{os.path.join(r,"movie.xyz")}}' filetype='xyz' 
                radius='2.0'  mfps='30' height=256 width=256
                style='border: 1px solid black; margin: 4px'></canvas>
			% end
		</p>

	</body>

</html>
