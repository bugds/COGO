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

		<title>Choosing referencial species</title>
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
				<input type="button" value="--&gt;" onclick="listbox_moveacross('notGood', 'good')">
				<br/>
				<input type="button" value="&lt;--" onclick="listbox_moveacross('good', 'notGood')">
			</td>
			<td>
				<p>Referencial countries</p>
				<form action="/calculate" method="post" onsubmit="placeInHidden(';', 'good', 'referencial');">
				    <input type="hidden" name="referencial" id="referencial">
				    <input type="hidden" name="eMail" value="{{eMail}}">
		            <input type="hidden" name="proteins" value="{{proteins}}">
					<select id="good" style="width: 200px;" size="20" multiple>
					</select>
					<input type="submit" value="Submit">
				</form>
			</td>
		</tr>
		</table>
	</body>
</html>