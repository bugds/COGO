<!DOCTYPE HTML>
<html>
	<head>
		<script type="text/javascript">
			function listbox_moveacross(sourceID, destID) {
				var src = document.getElementById(sourceID);
				var dest = document.getElementById(destID);

				for(var count=0; count < src.options.length; count++) {
					if(src.options[count].selected == true) {
							var option = src.options[count];

							var newOption = document.createElement("option");
							newOption.value = option.value;
							newOption.text = option.text;
							newOption.selected = true;
							try {
									 dest.add(newOption, null); //Standard
									 src.remove(count, null);
							 }catch(error) {
									 dest.add(newOption); // IE only
									 src.remove(count);
							 }
							count--;
					}
				}
			}

            function placeInHidden(delim, selStr, hidStr) {
                var selObj = document.getElementById(selStr);
                var hideObj = document.getElementById(hidStr);
                hideObj.value = '';

                for (var i=0; i<selObj.options.length; i++) {
                hideObj.value = hideObj.value ==
                    '' ? selObj.options[i].value : hideObj.value + delim + selObj.options[i].value;
                }
            }
		</script>

		<title>Listbox JavaScript functions</title>
	</head>

	<body>
		<table>
		<tr valign="top">
			<td>
				<p>Non-referencial countries</p>
				<select id="notGood" style="width: 200px;" size="20" multiple>{{!options}}
				</select>
			</td>
			<td align="center" valign="middle">
				<br/>
				<p id='toLeft' onclick="listbox_moveacross('notGood', 'good')">&gt;&gt;</p>
				<br/>
				<p id='toRight' onclick="listbox_moveacross('good', 'notGood')">&lt;&lt;</p>
			</td>
			<td>
				<p>Referencial countries</p>
				<form action="/calculate" method="get" onsubmit="placeInHidden(';', 'good', 'hidden');">
				    <input type="hidden" name="hidden" id="hidden">
					<select id="good" style="width: 200px;" size="20" multiple>
					</select>
					<input type="submit" value="Submit">
				</form>
			</td>
		</tr>
		</table>
	</body>
</html>