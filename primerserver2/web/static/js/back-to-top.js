$( document ).ready(function() {
			
			// appends html to body element - make sure you are linked to bootstrap icon fonts. 
			$( 'body' ).append('<a class="back-to-top" href="#" title="Back to Top" data-toggle="tooltip"><i class="glyphicon glyphicon-chevron-up"></i></a>');
  			$('[data-toggle="tooltip"]').tooltip({html: true});
			
		    var offset = 50; // sets distance in pixels when back to top will appear 
			var fadeintime = 400;  // fade in duration 
			$(window).scroll(function() {
				if ($(this).scrollTop() > offset) {
					$('.back-to-top').fadeIn(fadeintime);
				} else {
					$('.back-to-top').fadeOut(fadeintime);
				}
			});
			
			
			$('.back-to-top').click(function(){
			$("html, body").animate({scrollTop:0}, 'slow');	
				return false;
				
				});	

});	

