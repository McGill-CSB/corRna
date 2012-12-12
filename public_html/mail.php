<?php
ini_set('display_errors', '1');
	if(!empty($_POST['e'])) {
		if($_POST['e']!="") {
			$to = "kam.alfred@gmail.com";
			$subject = "[BUGS-Cor-RNA]";
			$message = $_POST["e"];
			mail($to, $subject, $message);
			echo "pass";
		}
	}
?>
