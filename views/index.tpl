<!DOCTYPE HTML>
<html>
<head>
<meta charset="utf-8">
<title>Input form</title>
</head>
<body>
    <form action="/reference" method="post" enctype="multipart/form-data">
        <p>Enter your e-mail:<input type="text" method="get" name="eMail"></p>
        <div>
            <label for="analysis">Choose file to submit...</label>
            <input type="file" id="analysis" name="analysis">
        </div>
        <div>
             <button>Submit</button>
        </div>
	<br>
        <div>
            <label for="reanalysis">...or upload saved analysis</label>
            <input type="file" id="reanalysis" name="reanalysis">
        </div>
        <div>
             <button>Upload</button>
        </div>
    </form>
</body>
</html>