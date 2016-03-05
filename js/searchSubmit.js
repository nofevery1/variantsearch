



$(document).bind("ajaxSend", function(){
   $("#loading").show();
 }).bind("ajaxComplete", function(){
   $("#loading").hide();
 });


$( "#variantData" ).submit(function( event ) {
  if ( $("#rsid").val().length === 0 && $( "#gene" ).val().length === 0 ) {
    $( ".required" ).addClass("has-error");
    event.preventDefault();
  }
  else if ( $("#rsid").val().length === 0 && $( "#proteinPosition").val().length ===0 ) {
    $( ".required" ).addClass("has-error");
    event.preventDefault();
  } else {
    console.log($("input#gene").val())
    console.log($("input#rsid").val())
    console.log($("input#proteinPosition").val())
    $( "#variantData" ).hide();
    $.ajax({
      url: "http://0.0.0.0:8080",
      data: {
        gene: $("input#gene").val(),
        rsid: $("input#rsid").val(),
        proteinPosition: $("input#proteinPosition").val()
      },
      type: "GET",
      dataType : "json",
      success: function( json ) {
        console.log("successs")
        console.log(json)
        var resultCount = $("<p>", { "class": "lead", "text":json.length+" papers found."});
        resultCount.appendTo("#searchForm");
        for(var i = 0; i < json.length; i++) {
          var obj = json[i];
          var paper = $( "<div>", { "class":"paper", "id": obj.pmid} );
          paper.appendTo("#searchForm");
          var address = $("<a>", { "html":"<strong>"+obj.title+"</strong>", "href":obj.address }).appendTo("#"+obj.pmid);
          var abstract = $("<p>", { "text":obj.abstract}).appendTo("#"+obj.pmid);
          console.log(obj.id);
        }
      },
      error: function( xhr, status, errorThrown ) {
        console.log( "Sorry, there was a problem!" );
        console.log( "Error: " + errorThrown );
        console.log( "Status: " + status );
        console.dir( xhr );
      },
      complete: function( xhr, status ) {
        console.log()
        console.log( "The request is complete!" );
      }
    });
    event.preventDefault();
  }
});
