<!DOCTYPE HTML>
<html>
	<head>
		<title>Choosing referencial genes</title>
	</head>

	<body>
        <form action="/calculate" method="post" target="_blank">{{!radioButtons}}
            <input type="hidden" name="eMail" value="{{eMail}}">
            <input type="hidden" name="proteins" value="{{proteins}}">
            <input type="hidden" name="blastDict" value="{{blastDict}}">
            <input type="hidden" name="referencial" value="{{referencial}}">
            <input type="submit" value="Submit">
        </form>
	</body>
</html>